mqtxt_to_eset <- function(filename,
						  normalize = TRUE,
						  impute = TRUE,
						  tidymods = TRUE,
						  remove_decoys = TRUE,
						  remove_contaminants = TRUE,
						  additional_filter = NULL,
						  transformfun = log2,
						  drop_cols = FALSE,
						  normalization_method = NULL,
						  missing_value = 0,
						  use_lfq = FALSE,
						  spread_site_multiplicity = FALSE) {
	# TODO consider changing this backend to fread
	esetData <- data.table::fread(filename, integer64 = 'double')

	if (remove_decoys) {
		if ("Reverse" %in% colnames(esetData)) {
			esetData <- dplyr::filter(
				esetData,
				is.na(esetData[["Reverse"]]) |
					esetData[["Reverse"]] == "" )
		} else {
			warning("No Reverse Column in the samples, skipping decoy removal\n")
		}
	}

	if (remove_contaminants) {
		if ("Potential contaminant" %in% colnames(esetData)) {
			esetData <- dplyr::filter(
				esetData,
				is.na(esetData[["Potential contaminant"]]) |
					esetData[["Potential contaminant"]] == "")
		} else {
			warning("No Potential contaminant Column in the samples, skipping decoy removal\n")
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


	if (spread_site_multiplicity) {
		abundance_cols <- esetData %>%
			dplyr::select(dplyr::matches("Intensity .+|^id$")) %>%
			dplyr::select(dplyr::matches("__|^id$")) %>%
			tidyr::gather(key, value, -id) %>%
			tidyr::extract(key, c("Sample", "Multiplicity"), "(I.+)___(.+)") %>%
			tidyr::spread(Sample, value)
		esetData <- left_join(
			abundance_cols,
			esetData %>%
				dplyr::select(-dplyr::matches("Intensity .+")))
		abundance_cols <- esetData %>%
			dplyr::select(dplyr::matches("Intensity .+"))
	} else {
		abundance_cols <- esetData %>%
			dplyr::select(dplyr::matches("Intensity.+")) %>%
			dplyr::select(-dplyr::matches("__"))
	}



	if (any(grepl("LFQ", colnames(abundance_cols)))) {
		if (use_lfq) {
			abundance_cols <- dplyr::select(dplyr::matches("LFQ"))
		} else {
			warning(
				paste0(
					"Removing columns with LFQ in the name, check if you",
					" want to use them with the 'use_lfq argument'"))
			abundance_cols <- dplyr::select(-dplyr::matches("LFQ"))
		}
	}

	abundance_cols[] <- lapply(abundance_cols, as.double)

	abundance_cols <- data.matrix(abundance_cols)


	if (!is.null(missing_value)) {
		# TEMPORARY FIX for cases where 64bit ints are converted to 10e-316
		abundance_cols[abundance_cols < 1] <- missing_value
		abundance_cols[abundance_cols == missing_value] <- NA
	}


	if (normalize) {
		abundance_cols <- limma::normalizeBetweenArrays(
			abundance_cols,
			method = normalization_method)
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
		# TODO add a way to add the min to each column by separate
		num_missing <- sum(is.na(Biobase::exprs(Eset)))
		warning(glue::glue(
				"Imputing {num_missing} Values from the dataset\n",
				num_missing))
		Biobase::exprs(Eset)[is.na(Biobase::exprs(Eset))] <-
			min(Biobase::exprs(Eset), na.rm = TRUE)
		Biobase::exprs(Eset)[0 > Biobase::exprs(Eset)] <-
			min(Biobase::exprs(Eset)[0 < Biobase::exprs(Eset)], na.rm = TRUE)
	}

	return(Eset)
}



#' Gets the Evidence entries for a series of ids
#'
#' @param evidence_ids Evidence Ids to be searched (can be a multi-numeric)
#' @param evidence_data data.frame, containing the evidence data (output of fread or equivalent)
#' @param multiplicity integer, vector indicating number of modifications of the type per peptide
#' @param evidence_multiplicity_col character, column name to use for the multiplicity in the evidence_data, defaults to "Phospho (STY)"
#' @param multi_numeric_separator character, separator for multi-numeric columns, defaults to ';'
#'
#' @return data.frame, containing the matching evidence entries
#' @export
#'
#' @examples
mq_site_to_evidence <- function(
	evidence_ids,
	evidence_data,
	multiplicity = NA,
	evidence_multiplicity_col = "Phospho (STY)",
	multi_numeric_separator = ';') {


	if (!all(is.na(multiplicity))) {
		site_glued_ids <- tibble::data_frame(
			ev_ids = strsplit(evidence_ids, ';'),
			Multiplicity = multiplicity) %>%
			tidyr::unnest() %>%
			dplyr::mutate(
				glued_ids = glue::glue(
					'{ev_ids}:{Multiplicity}',
					ev_ids = ev_ids,
					Multiplicity = Multiplicity )
			) %>%
			.[['glued_ids']]

		evidence_data <- evidence_data  %>%
			dplyr::mutate(
				glued_ids = glue::glue(
					'{ev_ids}:{Multiplicity}',
					ev_ids = evidence_data[['id']],
					Multiplicity = evidence_data[[evidence_multiplicity_col]] )
			)

		evidence_data <- evidence_data %>%
			dplyr::filter(glued_ids %in% site_glued_ids)
	} else {
		evidence_data <- evidence_data %>%
			dplyr::filter(
				id %in% unlist(strsplit(evidence_ids, ';'),
							   use.names = FALSE))
	}
	return(evidence_data)

}




#' Converts Maxquant to Skyline style peptide annotations
#'
#' @param annotated_peptides character, maxquant style annotated peptides \_PEPT(ph)IDE\_
#'
#' @return character, skyline style annotated peptide PEPT[+80]IDE
#' @export
#'
#' @examples
mq_to_skyline_peptides <- function(annotated_peptides) {
	# TODO add argument to convert other modifications to the skyline version
	warning("Only Phosphorylation, acetylation and oxidation have been implemented")
	annotated_sequences <- annotated_peptides %>%
		gsub('_', '', .) %>%
		strsplit(';') %>%
		gsub('\\(ph\\)', '[+80]', .) %>%
		gsub('\\(ox\\)', '[+16]', .) %>%
		gsub('\\(ac\\)', '[+42]', .) %>%
		unlist(use.names = FALSE)
	return(annotated_sequences)
}


mqtxt_to_sqlite <- function(txt_dir,
							db_name,
						  verbose = F,
						  light =  F) {
	stop("Not implemented yet")
}



make_exp_design <- function(
	eset,
	factor_vector = c(WT = 'WT', KO = 'KO'),
	reference_factor = 'KO',
	batch_vector = NULL) {

	get_factors.colnames <- function(my_colnames, named_vector) {
		# EXAMPLE
		# > mycolumns <- c('A1', 'A2', 'B1', 'B2')
		# > my_categories <- c('COOL_1' = 'A', 'COOL_2' = 'B')
		# > get_factors.colnames(mycolumns, my_categories)
		# [1] COOL_1 COOL_1 COOL_2 COOL_2
		# Levels: COOL_1 COOL_2
		matching_table <- sapply(named_vector, function(x){grepl(x, my_colnames)})

		if (any(rowSums(matching_table) != 1)) {
			print('Variables in the provided vector are not mutually exclusive\n')
			print(matching_table)
			stop('Variables in the provided vector are not mutually exclusive\n')
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

from_clipboard  <- function(header=TRUE, ...) {
	read.table("clipboard", sep = "\t",header = header, as.is = TRUE, ...)
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


count_missing_by.ExpressionSet <- function(data, named_character_vector) {
	mydata <- Biobase::exprs(data)
	ret <- count_missing_by.data.frame(
		mydata,
		named_character_vector = named_character_vector)
	return(ret)
}

count_missing_by.data.frame <- function(data, named_character_vector) {

	mynames <- colnames(data)
	column_selectors <- purrr::map(named_character_vector, function(x){ grepl(x, mynames) })

	df <- purrr::map_dfc(
		column_selectors,
		function(x){
			columns <- data[, x]
			if (is.null(dim(columns))) columns <- data.matrix(columns)
			rowSums(is.na(columns))})
	return(df)
}
