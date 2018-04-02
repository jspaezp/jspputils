library(testthat)
library(jspputils)

test_check("jspputils")


# setwd('I:/Elements (F drive)/Sebastian/2017/CDKL5 Project/09192017-Justine-InVitro')
#
# pdtxt_to_eset("012216-CDKL5-Justine_PeptideGroups.txt")
#
# myeset <- pdtxt_to_eset("012216-CDKL5-Justine_PeptideGroups.txt")
#
# myeset <- pdtxt_to_eset("012217-CDKL5-Justine_PeptideGroups.txt")
#
# design <- make_exp_design(myeset)
#
# fit <- lmFit(myeset, design)
# contr <- makeContrasts(WT-KO,levels=design)
# fit.contr <- eBayes(contrasts.fit(fit,contr))
#
# fit.contr
#
# maplot(fit.contr)
# foo <- topTable(fit.contr, number = Inf) %>%
# 	dplyr::filter( adj.P.Val < 0.05) %>%
# 	dplyr::filter(logFC > 2)
#
# foo %>%
# 	.$Master.Protein.Accessions %>%
# 	gsub("-[0-9]+", '', .) %>% grep("\\w+", . , value = TRUE) %>%
# 	writeLines("asdasdasdad.txt")
#
# volcanoplot(fit.contr)
#
# write.fit(fit.contr, file = 'writtenfit2.0.txt')
#
# ##### MQ
# #####
# myeset <- mqtxt_to_eset("I:/Elements (F drive)/Justine/CDKL5/In vitro mouse combined/combined/CDKL5 in vitro_LFQ_txt/18O-phospho (STY)Sites.txt", impute = TRUE)
#
# design <- make_exp_design(myeset)
#
# fit <- lmFit(myeset, design)
# contr <- makeContrasts(WT-KO,levels=design)
# fit.contr <- eBayes(contrasts.fit(fit,contr))
#
# fit.contr
# maplot(fit.contr)
# volcanoplot(fit.contr)
#
# dat <- read_tsv("I:/Elements (F drive)/Justine/CDKL5/In vitro mouse combined/combined/CDKL5 in vitro_LFQ_txt/18O-phospho (STY)Sites.txt")
#
# dat %>% dplyr::select(matches("Intensity CDKL5 KO")) %>% dplyr::select(matches("__"))  %>%
# 	(function(x){log10(x + 10000)}) %>% pairs()
#
# dat %>% dplyr::select(matches("Intensity CDKL5 WT")) %>% dplyr::select(matches("__"))  %>%
# 	(function(x){log10(x + 10000)}) %>% pairs()
#
# dat2 <- dat %>%
# 	dplyr::select(matches("Intensity CDKL5 WT|name")) %>%
# 	dplyr::select(matches("__1|name")) %>% as.data.frame()
# 	(function(x){log10(x + 10000)}) %>% pairs()
#
# ##### MQ comm vs IP
# #####
# myeset <- mqtxt_to_eset("I:/Elements (F drive)/modificationSpecificPeptides.txt", impute = TRUE)
#
# design <- make_exp_design(myeset)
#
# fit <- limma::lmFit(myeset, design)
# contr <- limma::makeContrasts(WT - KO,levels = design)
# fit.contr <- limma::eBayes(contrasts.fit(fit,contr))
#
# fit.contr
# maplot(fit.contr)
# volcanoplot(fit.contr)
#
# dat <- readr::read_tsv("I:/Elements (F drive)/modificationSpecificPeptides.txt")
# dat <-  dat %>% dplyr::filter(`18O-phospho (STY)` != 0)
# logdat <- dat %>%
# 	dplyr::mutate(`Int. CDKL5 kinase domain` = log10(`Intensity COMM` + 10000),
# 				  `Int. CDKL5 full length` = log10(`Intensity IP` + 10000))
#
# aest = ggplot2::aes(x = `Int. CDKL5 kinase domain`, y = `Int. CDKL5 full length` )
# logdat %>%
# 	ggplot2::ggplot(aest) +
# 	ggplot2::geom_point() +
# 	ggplot2::theme_bw() +
# 	ggrepel::geom_text_repel(
# 		logdat %>%
# 			dplyr::filter(`Int. CDKL5 kinase domain` > 4) %>%
# 			dplyr::filter(`Int. CDKL5 full length` > 4) %>%
# 			dplyr::filter(!grepl("Phos", Modifications)),
# 		mapping = aes(x = `Int. CDKL5 kinase domain`,
# 					  y =`Int. CDKL5 full length`,
# 					  label = Sequence),
# 		force = 4, alpha = 0.8)
#
# dat %>% dplyr::select(matches("Intensity CDKL5 WT|Sequence")) %>%
# 	dplyr::select(matches("__"))  %>% %>%
# 	(function(x){log10(x + 10000)}) %>% pairs()
#
# dat2 <- dat %>%
# 	dplyr::select(matches("Intensity CDKL5 WT|name")) %>%
# 	dplyr::select(matches("__1|name")) %>% as.data.frame()
# (function(x){log10(x + 10000)}) %>% pairs()
#
#
# # library(UniProt.ws)
# up <- UniProt.ws::UniProt.ws(taxId=9606)
#
# keys <- fit.contr$genes$`Master Protein Accessions`
# columns(up)
#
# columns <- c("ENTRY-NAME")
# kt <- "UNIPROTKB"
# res <- UniProt.ws::select(up, keys, columns, kt)
# res
#
