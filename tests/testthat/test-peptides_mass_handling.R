context("peptides_mass_handling")

test_that("Mass Diff Calcuations Work", {
    expect_equal(get_massdiff('S', 'T'), 14.01565)
	expect_equal(get_massdiff('T', 'S'), -14.01565)
	expect_equal(get_massdiff('L', 'I'), 0)
	expect_equal(get_massdiff('STLI', 'LITS', verbose = FALSE), 0)
	expect_output(get_massdiff('S', 'T', verbose = TRUE), "S --> T = 14.01565")
	expect_silent(get_massdiff('S', 'T', verbose = FALSE))
})


test_that("Mass Diffs are calculated correctly from CIGAR strings", {
    CIGAR_strings <- c(
    	'ST', '14', '1ST2IV', '1TN12', 'TS3QE3ST1TR1KQ',  '4QK3SE1TY1KR',
    	'3QH3ST1TL1KE', '3QIWF4TA1', 'SA4QV3ST1TV1')
    expect_equal(round(cigar_massdiff('4QK3SE1TY1KR', verbose = FALSE), digits = 4) , 132.0687)
    expect_equal(round(cigar_massdiff(CIGAR_strings, verbose = FALSE), digits = 4),
    			 c(14.0156,   0.0000,   0.0000,
    			   12.9953,  56.0011, 132.0687,
    			   36.0000, -83.9960, -32.9487))
out <-
"S --> T = 14.01565
NA --> NA = 0
S --> T = 14.01565
I --> V = -14.01565
T --> N = 12.99525
T --> S = -14.01565
Q --> E = 0.984009999999984
S --> T = 14.01565
T --> R = 55.05343
K --> Q = -0.0363799999999799
Q --> K = 0.0363799999999799
S --> E = 42.01056
T --> Y = 62.01565
K --> R = 28.00615
Q --> H = 9.00032999999999
S --> T = 14.01565
T --> L = 12.03638
K --> E = 0.947630000000004
Q --> I = -14.97452
W --> F = -39.0109
T --> A = -30.01057
S --> A = -15.99492
Q --> V = -28.99017
S --> T = 14.01565
T --> V = -1.97927"

    warn <- 'Either Base or replacement sequenes are NA'
    expect_warning({
    	expect_output({
    		round(cigar_massdiff(CIGAR_strings, verbose = TRUE), digits = 4)
    	}, out)}, warn)
})

test_cigar_massdiff <-  function() {
	cigar_massdiff('ST', verbose = FALSE)
	cigar_massdiff('TS', verbose = FALSE)
	cigar_massdiff('14', verbose = FALSE)
	cigar_massdiff('1ST2IV', verbose = FALSE)

	cigar_massdiff('1TN12', verbose = TRUE)
	cigar_massdiff('14', verbose = TRUE)
	cigar_massdiff('1ST2IV', verbose = TRUE)
	cigar_massdiff(c('1ST2IV', '14'), verbose = TRUE)
	cigar_massdiff(
		c('1TN12', 'TS3QE3ST1TR1KQ',
		  '4QK3SE1TY1KR', '3QH3ST1TL1KE',
		  '3QIWF4TA1', 'SA4QV3ST1TV1'),
		verbose = FALSE)
}
