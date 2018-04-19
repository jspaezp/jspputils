mqtxt_to_eset <- function(filename,
						  normalize = TRUE,
						  impute = TRUE,
						  tidymods = TRUE,
						  remove_decoys = TRUE,
						  remove_contaminants = TRUE,
						  additional_filter = NULL,
						  transformfun = log2,
						  drop_cols = FALSE,
						  missing_value = 0) {
	# TODO consider changing this backend to fread
	esetData <- readr::read_tsv(filename)

	if (remove_decoys) {
		if ("Reverse" %in% colnames(esetData)) {
			esetData <- dplyr::filter(
				esetData,
				is.na(esetData[["Reverse"]]))
		} else {
			warning("No Reverse Column in the samples, skipping decoy removal")
		}
	}

	if (remove_contaminants) {
		if ("Potential contaminant" %in% colnames(esetData)) {
			esetData <- dplyr::filter(
				esetData,
				is.na(esetData[["Potential contaminant"]]))
		} else {
			warning("No Potential contaminant Column in the samples, skipping decoy removal")
		}
	}

	if (is.character(additional_filter)) {
		esetData <- seplyr::filter_se(
			esetData,
			additional_filter)
	}

	if (is.character(drop_cols)) {
		esetData <- esetData[, !grepl(drop_cols, colnames(esetData))]
	}

	abundance_cols <- dplyr::select(esetData, dplyr::matches("Intensity.+")) %>%
		dplyr::select(-dplyr::matches("__")) %>%
		data.matrix()

	if (is.null(missing_value)) {
		abundance_cols[abundance_cols == missing_value] <- NA
	}


	if (normalize) {
		abundance_cols <- limma::normalizeBetweenArrays(abundance_cols)
	}

	abundance_cols <- abundance_cols %>% transformfun()
	abundance_cols[is.infinite(abundance_cols)] <- NA

	esetData$`Has Missing` <- is.na(rowSums(abundance_cols))
	esetData$`n.Missing` <- rowSums(is.na(abundance_cols))

	# Adds original intensities inside the feature data
	mycolnames <- colnames(esetData)
	feature_cols <- seplyr::rename_se(
		esetData,
		wrapr::named_map_builder(
			paste("Original", grep("Intensity", mycolnames, value = TRUE)),
			grep("Intensity", mycolnames, value = TRUE)))

	Eset <- Biobase::ExpressionSet(
		abundance_cols,
		featureData = as(feature_cols, "AnnotatedDataFrame"))

	if (impute) {
		warning("Imputing Values from the dataset")
		Biobase::exprs(Eset)[is.na(Biobase::exprs(Eset))] <-
			min(Biobase::exprs(Eset), na.rm = TRUE)
		Biobase::exprs(Eset)[0 > Biobase::exprs(Eset)] <-
			min(Biobase::exprs(Eset)[0 < Biobase::exprs(Eset)], na.rm = TRUE)
	}

	# TODO add warning to say what was imputed

	return(Eset)
}



### Design Maker
#' Title
#'
#' @param eset
#' @param factor_vector
#'
#' @return
#' @export
#'
#' @examples
make_exp_design <- function(eset, factor_vector = c(WT = 'WT', KO = 'KO'), reference_factor = 'KO', batch_vector = NULL) {
	# Makes only simple single factor designs
	# List should be in the form of c(WT = 'WT', KO = 'KO'), being WT and KO
	# part of the names of the columns in the given eset
	#

	get_factors.colnames <- function(my_colnames, named_vector) {
		# EXAMPLE
		# > mycolumns <- c('A1', 'A2', 'B1', 'B2')
		# > my_categories <- c('COOL_1' = 'A', 'COOL_2' = 'B')
		# > get_factors.colnames(mycolumns, my_categories)
		# [1] COOL_1 COOL_1 COOL_2 COOL_2
		# Levels: COOL_1 COOL_2
		matching_table <- sapply(named_vector, function(x){grepl(x, my_colnames)})

		if (any(rowSums(matching_table) != 1)) {
			print(matching_table)
			stop('Variables in the provided vctor are not mutually exclusive')
		}

		s_mainfactor <- factor(levels = colnames(matching_table))

		for (i in colnames(matching_table)) {
			s_mainfactor[(matching_table[,i])] <- i
		}
		return(s_mainfactor)
	}

	s_mainfactor <- get_factors.colnames(colnames(eset), factor_vector)


	if (is.character(reference_factor)) {
		print(paste("Using", reference_factor, "as a reference"))
		s_mainfactor <- relevel(s_mainfactor, reference_factor)
	}

	if (is.null(batch_vector)) {
		design <- model.matrix(~0+s_mainfactor)
		colnames(design) <- levels(s_mainfactor)
	} else {
		s_batch <- get_factors.colnames(colnames(eset), batch_vector)
		design <- model.matrix(~ s_batch + s_mainfactor)
	}

	rownames(design) <- colnames(eset)
	return(design)
}

# TODO make function to summarize the datasets ...
# number of samples,
# num of peptides/proteins/psm/msms/modified vs unmodified per sample
#

impute_eset <- function(eset, nPcs, seed = 6, ...){
	#myeset2 <- myeset[rowSums(data.matrix(Biobase::exprs(myeset)), na.rm = TRUE) != 0,]
	#Biobase::exprs(myeset2) <- Biobase::exprs(myeset2) - min(Biobase::exprs(myeset2), na.rm = TRUE)

	zero_variance_rows <- rowSums(
		data.matrix(Biobase::exprs(eset)), na.rm = TRUE
		) == 0
	if (sum(zero_variance_rows) > 0) {
		warning(paste0("Found rows without any value (or all 0), they will be removed -",
					   sum(zero_variance_rows),
					   "- rows"))
		eset <- eset[!zero_variance_rows,]
}

	pc <- pcaMethods::pca(eset, nPcs = nPcs, method = "ppca", seed = seed, ...)
	#pc <- pca(myeset2, nPcs = 3,  method = 'bpca')
	imputed <- t(pcaMethods::completeObs(pc))
	#pairs(imputed)
	Biobase::exprs(eset) <- imputed
	return(eset)
}


to_clipboard <- function(x) {
	write.table(x, "clipboard-512", sep = "\t", row.names = FALSE)
	print("Done")
}

read_panther <- function(file, verbose = TRUE) {

	raw_lines <- readLines(file)

	tablelines <- raw_lines %>% grepl("^.*(\\t).*(\\1).*$", .)
	contextlines <- raw_lines[{raw_lines %>% grepl("^.+$", .)} & !tablelines]

	context <- strsplit(contextlines, '\t')

	pantherdf <- readr::read_tsv(
		file,
		skip =  min(which(tablelines)) - 1)

	foldcolumn <- colnames(pantherdf) %>%
		grep("fold",
			 .,
			 ignore.case = TRUE,
			 value = TRUE)

	pantherdf[[foldcolumn]] <- gsub("[<> ]", "", pantherdf[[foldcolumn]]) %>%
		as.numeric()

	for (i in context) {
		attr(pantherdf, make.names(i[[1]])) <- i[[2]]
	}

	if (verbose) {
		print(attributes(pantherdf))
	}

	return(pantherdf)
}


count_missing_by.data.frame <- function(data, named_character_vector) {

	mynames <- colnames(data)
	column_selectors <- purrr::map(named_character_vector, function(x){ grepl(x, mynames) })

	df <- purrr::map_dfc(
		column_selectors,
		function(x){ rowSums(is.na(data[, x])) })
	return(df)
}
