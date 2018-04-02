mqtxt_to_eset <- function(filename,
						  normalize = TRUE,
						  impute = TRUE,
						  tidymods = TRUE,
						  remove_decoys = TRUE,
						  remove_contaminants = TRUE,
						  additional_filter = NULL,
						  transformfun = log2) {
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

	abundance_cols <- dplyr::select(esetData, dplyr::matches("Intensity.+")) %>%
		dplyr::select(-dplyr::matches("__")) %>%
		data.matrix() %>%
		transformfun()

	abundance_cols[is.infinite(abundance_cols)] <- NA
	esetData$`Has Missing` <- is.na(rowSums(abundance_cols))
	esetData$`n.Missing` <- rowSums(is.na(abundance_cols))

	feature_cols <- dplyr::select(esetData, -dplyr::matches("Intensity"))

	Eset <- Biobase::ExpressionSet(
		abundance_cols,
		featureData = as(feature_cols, "AnnotatedDataFrame"))

	if (impute) {
		Biobase::exprs(Eset)[is.na(Biobase::exprs(Eset))] <-
			min(Biobase::exprs(Eset), na.rm = TRUE)
		Biobase::exprs(Eset)[0 > Biobase::exprs(Eset)] <-
			min(Biobase::exprs(Eset)[0 < Biobase::exprs(Eset)], na.rm = TRUE)
	}

	if (normalize) {
		Biobase::exprs(Eset) <- limma::normalizeBetweenArrays(Biobase::exprs(Eset))
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
make_exp_design <- function(eset, factor_vector = c(WT = 'WT', KO = 'KO')) {
	# Makes only simple single factor designs
	# List should be in the form of c(WT = 'WT', KO = 'KO'), being WT and KO
	# part of the names of the columns in the given eset

	joint_factors <- paste(factor_vector, sep = '|', collapse = '|')
	cleaner_regex <- paste0(".*(", joint_factors, ").*$")
	samples <- gsub(cleaner_regex, "\\1", colnames(eset))

	s <- factor(factor_vector[samples])
	design <- model.matrix(~0+s)
	colnames(design) <- levels(s)
	return(design)
}

# TODO make function to summarize the datasets ...
# number of samples,
# num of peptides/proteins/psm/msms/modified vs unmodified per sample
#

impute_eset <- function(eset, nPcs){
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

	pc <- pcaMethods::pca(eset, nPcs = nPcs, method = "ppca")
	#pc <- pca(myeset2, nPcs = 3,  method = 'bpca')
	imputed <- t(pcaMethods::completeObs(pc))
	#pairs(imputed)
	Biobase::exprs(eset) <- imputed
	return(eset)
}


plotpca <- function(eset, pcaresult) {
	df <- merge(pcaMethods::scores(pcaresult), Biobase::pData(eset), by = 0)
	ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = as.character(Row.names))) +
		ggplot2::geom_point() +
		ggplot2::xlab(paste("PC1", pcaresult@R2[1] * 100, "% of variance")) +
		ggplot2::ylab(paste("PC2", pcaresult@R2[2] * 100, "% of variance"))
}
