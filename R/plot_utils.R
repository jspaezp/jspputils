
ggplot <- function(..., xlab_name, ylab_name) {
	g <- ggplot2::ggplot(...) +
		ggplot2::theme_bw() +
		ggplot2::labs(x = xlab_name, y = ylab_name)
	return(g)
}

# TODO chech for duplication in maplot
pvaluecolours <- function(name = "P. Value", ...){
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
			title = name,
			draw.ulim = FALSE,
			nbin = 22),
		...)
}


extract_fit_table <- function(fit) {
	fittable <- limma::topTable(fit, number = Inf)
	fittable$ID <- rownames(fittable)
	colnames(fittable) <- gsub("^genes\\.", "", colnames(fittable))

	fullfitdf <- as.data.frame(fit)
	fullfitdf$ID <- rownames(fullfitdf)
	colnames(fullfitdf) <- gsub("^genes\\.", "", colnames(fullfitdf))

	#nonnumeric_cols <- colnames(fittable)[sapply(fittable, function(x){!is.numeric(x)})]
	#joincols <- intersect(intersect(colnames(fullfitdf) , colnames(fittable) ), nonnumeric_cols)

	fullfitdf2 <- dplyr::right_join(
		fullfitdf,
		fittable
		) # by  = joincols

	if (nrow(fullfitdf2) != nrow(fittable) | nrow(fullfitdf2) != nrow(fullfitdf)) {
		messagetemplate <- c("Number of rows (",
							 ") in the original and the generated table do not match,",
							 " try manually extracting the table")
		print(c(nrow(fittable), nrow(fullfitdf)))
		message <- paste(messagetemplate[[1]], nrow(fittable), nrow(fullfitdf), messagetemplate[[2]])
		stop(message)
	}
	return(fullfitdf2)
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

volcanoplot.data.frame <-  function(fullfitdf, foldchange = "logFC",
									pval = "adj.P.Val", nlabel = 5,
									labelcol = "Gene.Names",
									add_aes_base = NULL,
									add_aes_label = add_aes_base,
									ylab_name = "PLEASE LABEL YOUR AXES",
									xlab_name = "PLEASE LABEL YOUR AXES",
									...) {
	# TODO add something to extract the name of the contrasts and append it to
	# the logFC name

	aes_list <- list(
		x = foldchange,
		y = pval,
		label = labelcol)

	tmp_aes <- dispatch_aes(aes_list, add_aes_base)
	tmp_aes_label <- dispatch_aes(aes_list, add_aes_label)

	g <- ggplot(fullfitdf, tmp_aes, xlab_name = xlab_name, ylab_name = ylab_name) +
		ggplot2::scale_y_log10(breaks = c(1, 0.5,0.1, 0.05,0.01,0.05 , 0.01)) +
		ggplot2::geom_point(alpha = 0.3)
	g <- g + labeltopn(fullfitdf, n = nlabel,
					   mapping = tmp_aes_label, arrangeTerms = pval, ...)
	return(g)

}

# TODO add a way to accept a multi-contrast fit as an argument
volcanoplot <- function(fit, ...) {
	fullfitdf <- extract_fit_table(fit)
	g <- volcanoplot.data.frame(fullfitdf, ...)
	return(g)
}

maplot.data.frame <- function(fullfitdf, meanexpression = 'Amean', foldchange = "logFC",
							  pval = "adj.P.Val", showmissing = TRUE,
							  nlabel = 10,
							  colour = pval,
							  colour_limits = c(0,0.1),
							  add_aes_base = NULL,
							  add_aes_label = add_aes_base,
							  labelcol = "Gene.Names",
							  ylab_name = "PLEASE LABEL YOUR AXES",
							  xlab_name = "PLEASE LABEL YOUR AXES",
							  label_arrangeTerms = paste0("desc(abs(", foldchange,"))"),
							  ...) {
	if (showmissing) {
		maxfoldchange <- 1 + max(fullfitdf[[foldchange]], na.rm = TRUE)
		fullfitdf[[foldchange]][is.na(fullfitdf[[foldchange]])] <- maxfoldchange
	}

	aes_list <- list(x = meanexpression, y = foldchange,
					 colour = colour, label = labelcol )

	tmp_aes <- dispatch_aes(aes_list, add_aes_base)
	tmp_aes_label <- dispatch_aes(aes_list, add_aes_label)


	g <- ggplot(fullfitdf,
				tmp_aes,
				xlab_name = xlab_name,
				ylab_name = ylab_name) +
		ggplot2::geom_point(alpha = 0.8) +
		ggplot2::geom_point(
			data = subset(
				fullfitdf,
				colour_limits[1] < fullfitdf[[colour]] &
					colour_limits[2] > fullfitdf[[colour]])) +
		ggplot2::scale_colour_gradientn(
			colours = c("red",
						"blue"),
			# values = c(0, 0.1),
			# labels = c(0.001, 0.05, 1),
			# breaks = c(0, 0.1),
			# trans = "sqrt",
			space = "Lab",
			na.value = "gray",
			limits = colour_limits,
			guide = ggplot2::guide_colourbar(
				title = "P. Value",
				draw.ulim = FALSE,
				nbin = 22))

	ifelse("filterterms" %in% names(list(...)),
		   no = {filterterms <- paste0(pval, " < 0.05 ") },
		   yes = {filterterms <- list(...)[["filterterms"]] })
	g <- g + labeltopn(fullfitdf, n = nlabel,
					   mapping = tmp_aes_label,
					   arrangeTerms = label_arrangeTerms,
					   filterterms =  filterterms)

	return(g)
}

maplot <- function(fit,  ...) {
	fullfitdf <- extract_fit_table(fit)
	g <- maplot.data.frame(fullfitdf, ...)
	return(g)
}

plot_ranked_folds.data.frame <- function(fullfitdf, nlabel = 10,
										 folds = "logFC",
										 labelcol = "Gene.Names",
										 add_aes_base = NULL,
										 add_aes_label = add_aes_base,
										 ylab_name = "PLEASE LABEL YOUR AXES",
										 xlab_name = "PLEASE LABEL YOUR AXES",
										 ...) {
	fullfitdf[["rank"]] <- rank(fullfitdf[[folds]], ties.method = "random")

	aes_list <- list( y = folds, x = "rank", label = labelcol, ...)

	tmp_aes <- dispatch_aes(aes_list, add_aes_base)
	tmp_aes_label <- dispatch_aes(aes_list, add_aes_label)

	g <- ggplot(fullfitdf, tmp_aes, xlab_name = xlab_name, ylab_name = ylab_name) +
		ggplot2::geom_point()
	g <- g + labeltopn(fullfitdf, n = nlabel,
					   mapping = tmp_aes_label,
					   arrangeTerms = paste0("desc(abs(", folds,"))"))
	return(g)
}

plot_ranked_folds <- function(fit, ...) {
	fullfitdf <- extract_fit_table(fit)
	g <- plot_ranked_folds.data.frame(fullfitdf, ...)
	return(g)

}


plot_dists.eset <- function(eset){
	Biobase::exprs(eset) %>%
		reshape2::melt() %>%
		ggplot(ggplot2::aes_string(x = "value", fill = "Var2"),
			   ylab_name = "Frequency", xlab_name = "value") +
		ggplot2::geom_density(alpha = 0.3)
}


plotpca <- function(eset, pcaresult) {
	df <- merge(pcaMethods::scores(pcaresult), Biobase::pData(eset), by = 0)
	ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = as.character(Row.names))) +
		ggplot2::geom_point() +
		ggplot2::xlab(paste("PC1", pcaresult@R2[1] * 100, "% of variance")) +
		ggplot2::ylab(paste("PC2", pcaresult@R2[2] * 100, "% of variance"))
}


plot_n_saveplotly <- function(ggplotobject, filename, ...) {
	plotlyobj <- plotly::ggplotly(ggplotobject)
	htmlwidgets::saveWidget(plotlyobj, file = filename, ...)
	print(ggplotobject)
	return(TRUE)
}

makealltheplots <- function(fit,
							coef = 's_mainfactorWT',
							contrast_name = coef,
							plotprefix = format(Sys.time(), "%Y%m%d_%H%M%S_"),
							nlabels = 20) {
	mydf <- limma::topTable(
		fit, coef = coef,
		number = Inf,
		confint = TRUE)

	mydf[[names_column]] <- paste0(
		mydf[[names_column]], ": ",
		mydf[["Positions.within.proteins"]])


	mydf[["Site.Names"]] <- paste0(
		mydf[[names_column]], "-",
		mydf[['Protein.names']],
		": ",
		mydf[["Positions.within.proteins"]])

	g <- maplot.data.frame(
		mydf,
		pval = "P.Value",
		foldchange = 'logFC',
		nlabel = nlabels,
		labelcol = "Site.Names",
		colour_limits = c(0,0.05),
		xlab_name = 'Mean log2 Intensity',
		ylab_name = paste0('Log2 Fold Change', contrast_name),
		add_aes_base = list(x = 'AveExpr'))
	plot_n_saveplotly(g, paste(plotprefix, 'maplot_longnames.html'))

	g <- maplot.data.frame(
		mydf,
		pval = "P.Value",
		foldchange = 'logFC',
		nlabel = nlabels,
		labelcol = names_column,
		colour_limits = c(0,0.05),
		add_aes_base = list(x = 'AveExpr'),
		xlab_name = 'Mean log2 Intensity',
		ylab_name = paste0('Log2 Fold Change', contrast_name),
		label_arrangeTerms = 'P.Value')
	plot_n_saveplotly(g, paste(plotprefix, 'maplot_shortnames.html'))

	g <- volcanoplot.data.frame(
		mydf,
		pval = 'P.Value',
		add_aes_base = list(colour = "as.character(n.Missing)"),
		nlabel = nlabels,
		labelcol = names_column,
		xlab_name = paste0('Log2 Fold Change', contrast_name),
		ylab_name = 'P.Value') +
		ggplot2::guides(colour = ggplot2::guide_legend(title =  "Number of \nMissing Values"))
	plot_n_saveplotly(g, paste(plotprefix, 'volcano_missingcolours.html'))

	g <- plot_ranked_folds.data.frame(
		mydf,
		labelcol = names_column,
		nlabel = nlabels,
		add_aes_base = list(colour = "as.character(n.Missing)"),
		xlab_name =  paste0('Ranked Log Fold Change', contrast_name),
		ylab_name =  paste0('Log2 Fold Change', contrast_name)) +
		ggplot2::guides(colour = ggplot2::guide_legend(title =  "Number of \nMissing Values"))
	plot_n_saveplotly(g, paste(plotprefix, 'rankedfolds_missingcolours.html'))


	g <- plot_ranked_folds.data.frame(
		mydf,
		labelcol = names_column,
		nlabel = nlabels,
		xlab_name = paste0('Ranked Log2 Fold Change', contrast_name),
		ylab_name = paste0('Log2 Fold Change Change', contrast_name))
	plot_n_saveplotly(g, paste(plotprefix, 'rankedfolds_nocolor.html'))


	g <- volcanoplot.data.frame(
		mydf,
		pval = 'P.Value',
		labelcol = names_column,
		nlabel = nlabels,
		xlab_name = paste0('Log2 Fold Change', contrast_name),
		ylab_name = 'P. Value')
	plot_n_saveplotly(g, paste(plotprefix, 'volcano_nocolor.html'))
	return(TRUE)
}
