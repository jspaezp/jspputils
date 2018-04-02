
extract_fit_table <- function(fit) {
	fittable <- limma::topTable(fit, number = Inf)
	fittable$ID <- rownames(fittable)

	fullfitdf <- as.data.frame(fit)
	fullfitdf$ID <- rownames(fullfitdf)

	fullfitdf <- dplyr::full_join(fullfitdf, fittable, by = c('ID' = 'ID'))
	return(fullfitdf)
}


labeltopn <- function(df, n=10, mapping, arrangeTerms, filterterms = NULL, ...) {
	df <- seplyr::filter_se(df, filterterms)
	newdf <- seplyr::arrange_se(df, arrangeTerms = arrangeTerms) %>%
		head(n)
	# print(newdf)
	labs <- ggrepel::geom_label_repel(
		data = newdf,
		mapping = mapping,
		min.segment.length = 0)
	return(labs)
}


volcanoplot <- function(fit, foldchange = "logFC",
						pval = "adj.P.Val", nlabel = 5,
						labelcol = "Gene.Names", ...) {

	fullfitdf <- extract_fit_table(fit)

	# TODO add something to extract the name of the contrasts and append it to
	# the logFC name
	tmp_aes <- ggplot2::aes_string(
		x = foldchange,
		y = pval,
		label = labelcol,
		...)
	g <- ggplot2::ggplot(fullfitdf, tmp_aes) +
		ggplot2::scale_y_log10(breaks = c(1, 0.5,0.1, 0.05,0.01,0.05 , 0.01)) +
		ggplot2::geom_point(alpha = 0.3) +
		ggplot2::theme_bw()

	if (nlabel > 0) {
		g <- g + labeltopn(fullfitdf, n = nlabel,
						   mapping = tmp_aes, arrangeTerms = pval)
	}
	return(g)

}


maplot <- function(fit, foldchange = "logFC",
				   pval = "adj.P.Val", showmissing = TRUE,
				   nlabel = 10, labelcol = "Gene.Names", ...) {
	fullfitdf <- extract_fit_table(fit)

	# if (showmissing) {
	# 	maxfoldchange <- 1 + max(fullfitdf[[foldchange]], na.rm = TRUE)
	# 	fullfitdf[]
	# }

	tmp_aes <- ggplot2::aes_string(x = "Amean",
								   y = foldchange,
								   colour = pval,
								   label = labelcol )

	g <- ggplot2::ggplot(fullfitdf,
						 tmp_aes) +
		ggplot2::geom_point(alpha = 0.8) +
		ggplot2::theme_bw() +
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

	if (nlabel > 0) {
		g <- g + labeltopn(fullfitdf, n = nlabel,
						   mapping = tmp_aes,
						   arrangeTerms = paste0("desc(abs(", foldchange,"))"),
						   filterterms = paste0(pval, " < 0.05"))
	}
	return(g)
}


plot_ranked_folds <- function(fit, nlabel = 10,
							  folds = "logFC", labelcol = "Gene.Names") {
	fullfitdf <- extract_fit_table(fit)
	fullfitdf[["rank"]] <- rank(fullfitdf[[folds]])

	tmp_aes <- ggplot2::aes_string(
		y = folds,
		x = "rank",
		label = labelcol)

	g <- ggplot2::ggplot(fullfitdf, tmp_aes) +
		ggplot2::geom_point() +
		ggplot2::theme_bw()
	if (nlabel > 0) {
		g <- g + labeltopn(fullfitdf, n = nlabel,
						   mapping = tmp_aes,
						   arrangeTerms = paste0("desc(abs(", folds,"))"))
	}
	return(g)
}


plot_dists.eset <- function(eset){
	Biobase::exprs(eset) %>%
		reshape2::melt() %>%
		ggplot2::ggplot(ggplot2::aes_string(x = "value", fill = "Var2")) +
		ggplot2::geom_density(alpha = 0.3) +
		ggplot2::theme_bw()
}
