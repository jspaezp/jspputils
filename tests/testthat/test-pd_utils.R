context("pd_utils")

test_that("multiplication works", {
	testingset <- c(
	"1xCarbamidomethyl [C3]; 1xPhospho [S15]; 1xPhospho18O [S14]",
	"1xCarbamidomethyl [C3]; 2xPhospho [S11; S15]; 1xPhospho18O [S9]",
	"1xPhospho [T16]; 1xPhospho18O [T29]",
	"1xPhospho18O [S]",
	"1xPhospho18O [S3(90.1)]",
	"3xCarbamidomethyl [C7; C14; C33]; 1xPhospho [S15]; 2xPhospho18O [S13; S19]",
	"1xPhospho [S/T]; 1xPhospho18O [S/T]",
	"1xPhospho [S12(97.1)]; 1xPhospho18O [S9(81.1)]",
	"1xOxidation [M23]; 2xPhospho18O [T9; S21]",
	"1xPhospho18O [Y8(100)]",
	""
	)
	split_mods(testingset)
	x <- get_mods(testingset)
	x
	purrr::map(.x = x, purrr::partial(find_mod, x = split_mods(testingset)[[1]]))
	mapmods(x, split_mods(testingset)[[1]])
	tidyfy_mods(testingset)
})
