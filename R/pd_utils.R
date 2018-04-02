
### Modition cleaning

split_mods <- function(x) {
	x <- strsplit(x, split = "(?<=.)(?=( [0-9]+x))", perl = TRUE)
	cleanws <- function(x) gsub('^\\s+|\\s+$|;$', '', x = x, perl = TRUE)
	x <- purrr::map(x, cleanws)
	return(x)
}

get_mods <- function(x) {
	splitmods <- unlist(split_mods(x), use.names = FALSE)
	cleanmods <- gsub('^\\s*[0-9]+x(\\w+) \\[.*$', '\\1', x, perl = TRUE)

	return(unique(cleanmods))
}


find_mod <- function(mod, x) {
	x <- grep(pattern = paste0(mod, ' '), x = x, value = TRUE)
	x <- ifelse(length(x) > 0, yes = x, no = NA)
	return(x)
}


mapmods <- function(mods, x){
	x <- purrr::map(.x = mods, purrr::partial(find_mod, x = x))
	names(x) <- make.names(mods)
	return(x)

}


tidyfy_mods <- function(x, keep_names = FALSE) {
	uniquemods <- get_mods(x)
	splitmods <- split_mods(x)

	df <- purrr::map_df(splitmods, purrr::partial(mapmods, mods = uniquemods))

	if (!keep_names) {
		df <- lapply(df, function(x) {gsub('^.*\\[(.*)\\]', '\\1', x)} )
		df <- data.frame(df)
	}

	return(df)
}

### Parser from PD exported txt file to an annotation object
pdtxt_to_eset <- function(filename, normalize = TRUE, impute = TRUE, tidymods = TRUE) {
	esetData <- readr::read_tsv(filename)

	abundance_cols <- dplyr::select(esetData, dplyr::matches("Abundance")) %>%
		data.matrix()

	feature_cols <- dplyr::select(esetData, -dplyr::matches("Abundance"))

	if (tidymods) {
		feature_cols <- cbind(
			tidyfy_mods(feature_cols$Modifications),feature_cols
		)
	}

	Eset <- Biobase::ExpressionSet(
		abundance_cols,
		featureData = as(feature_cols, "AnnotatedDataFrame"))

	Biobase::exprs(Eset) <- log(Biobase::exprs(Eset))

	if (normalize) {
		# TODO add a way to change the normalization method
		Biobase::exprs(Eset) <- dplyr::normalizeBetweenArrays(Biobase::exprs(Eset))
	}

	if (impute) {
		Biobase::exprs(Eset)[is.na(Biobase::exprs(Eset))] <-
			min(Biobase::exprs(Eset), na.rm = TRUE)
	}

	return(Eset)
}
