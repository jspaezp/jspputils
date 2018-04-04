
ggplot <- function(...) {
	g <- ggplot2::ggplot(...) +
		ggplot2::theme_bw()
	return(g)
}


extract_fit_table <- function(fit) {
	fittable <- limma::topTable(fit, number = Inf)
	fittable$ID <- rownames(fittable)
	colnames(fittable) <- gsub("^genes\\.", "", colnames(fittable))

	fullfitdf <- as.data.frame(fit)
	fullfitdf$ID <- rownames(fullfitdf)
	colnames(fullfitdf) <- gsub("^genes\\.", "", colnames(fullfitdf))

	fullfitdf <- dplyr::full_join(fullfitdf, fittable)
	if (nrow(fittable) != nrow(fullfitdf)) {
		stop("Number of rows does not match, try manually extracting the table")
	}
	return(fullfitdf)
}


labeltopn <- function(df, n=10, mapping, arrangeTerms, filterterms = NULL, ...) {
	if (n < 1) {
		return(NULL)
	}
	df <- seplyr::filter_se(df, filterterms)
	newdf <- seplyr::arrange_se(df, arrangeTerms = arrangeTerms) %>%
		head(n)
	# print(newdf)
	labs <- ggrepel::geom_label_repel(
		data = newdf,
		mapping = mapping,
		min.segment.length = 0, ...)
	return(labs)
}


dispatch_aes <- function(base_aes_list, add_aes_list = NULL) {
	base_aes_list <- c(base_aes_list, add_aes_list)
	base_aes_list <- rev(base_aes_list)[unique(names(base_aes_list))]
	aes <- do.call(ggplot2::aes_string, base_aes_list)

	return(aes)
}


volcanoplot <- function(fit, foldchange = "logFC",
						pval = "adj.P.Val", nlabel = 5,
						labelcol = "Gene.Names",
						add_aes_base = NULL,
						add_aes_label = add_aes_base,
						...) {

	fullfitdf <- extract_fit_table(fit)

	# TODO add something to extract the name of the contrasts and append it to
	# the logFC name

	aes_list <- list(
		x = foldchange,
		y = pval,
		label = labelcol)

	tmp_aes <- dispatch_aes(aes_list, add_aes_base)
	tmp_aes_label <- dispatch_aes(aes_list, add_aes_label)


	g <- ggplot(fullfitdf, tmp_aes) +
		ggplot2::scale_y_log10(breaks = c(1, 0.5,0.1, 0.05,0.01,0.05 , 0.01)) +
		ggplot2::geom_point(alpha = 0.3)
	g <- g + labeltopn(fullfitdf, n = nlabel,
					   mapping = tmp_aes_label, arrangeTerms = pval, ...)
	return(g)

}


maplot <- function(fit, foldchange = "logFC",
				   pval = "adj.P.Val", showmissing = TRUE,
				   nlabel = 10, colour = pval,
				   add_aes_base = NULL,
				   add_aes_label = add_aes_base,
				   labelcol = "Gene.Names", ...) {
	fullfitdf <- extract_fit_table(fit)

	if (showmissing) {
		maxfoldchange <- 1 + max(fullfitdf[[foldchange]], na.rm = TRUE)
		fullfitdf[[foldchange]][is.na(fullfitdf[[foldchange]])] <- maxfoldchange
	}

	aes_list <- list(x = "Amean", y = foldchange,
					 colour = colour, label = labelcol )

	tmp_aes <- dispatch_aes(aes_list, add_aes_base)
	tmp_aes_label <- dispatch_aes(aes_list, add_aes_label)


	g <- ggplot(fullfitdf,
				tmp_aes) +
		ggplot2::geom_point(alpha = 0.8) +
		ggplot2::scale_colour_gradientn(
			colours = c("red",
						"blue"),
			# values = c(0, 0.1),
			# labels = c(0.001, 0.05, 1),
			# breaks = c(0, 0.1),
			# trans = "sqrt",
			space = "Lab",
			na.value = "gray",
			limits = c(0,0.1),
			guide = ggplot2::guide_colourbar(
				title = "Adjusted\nP. Value\n",
				draw.ulim = FALSE,
				nbin = 22)) +
		ggplot2::labs(x = "Mean log Intensity", y = "Log Fold Change")

	ifelse("filterterms" %in% names(list(...)),
		   no = {filterterms <- paste0(pval, " < 0.05 ") },
		   yes = {filterterms <- list(...)[["filterterms"]] })
	g <- g + labeltopn(fullfitdf, n = nlabel,
					   mapping = tmp_aes_label,
					   arrangeTerms = paste0("desc(abs(", foldchange,"))"),
					   filterterms =  filterterms)

	return(g)
}


plot_ranked_folds <- function(fit, nlabel = 10,
							  folds = "logFC",
							  labelcol = "Gene.Names",
							  add_aes_base = NULL,
							  add_aes_label = add_aes_base,
							  ...) {
	fullfitdf <- extract_fit_table(fit)
	fullfitdf[["rank"]] <- rank(fullfitdf[[folds]], ties.method = "random")

	aes_list <- list( y = folds, x = "rank", label = labelcol, ...)

	tmp_aes <- dispatch_aes(aes_list, add_aes_base)
	tmp_aes_label <- dispatch_aes(aes_list, add_aes_label)

	g <- ggplot(fullfitdf, tmp_aes) +
		ggplot2::geom_point()
	g <- g + labeltopn(fullfitdf, n = nlabel,
					   mapping = tmp_aes_label,
					   arrangeTerms = paste0("desc(abs(", folds,"))"))
	return(g)
}


plot_dists.eset <- function(eset){
	Biobase::exprs(eset) %>%
		reshape2::melt() %>%
		ggplot(ggplot2::aes_string(x = "value", fill = "Var2")) +
		ggplot2::geom_density(alpha = 0.3)
}
