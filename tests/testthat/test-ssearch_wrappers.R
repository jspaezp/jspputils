context("ssearch_wrappers")

test_that("Parsing ssearch works", {

	expect_s3_class({
		sample_out <- c(
			"# C:\\Users\\SEBAST~1\\opt\\fasta36\\SSEARC~1.EXE -m 8C -E 0.01 C:\\Users\\SEBAST~1\\AppData\\Local\\Temp\\RtmpGYHOVK\\2b019357d09.fasta C:\\Users\\Sebastian\\Downloads\\uniprot-taxonomy-lambda.fasta",
			"# SSEARCH 36.3.8g Dec, 2017", "# Query: AAAAAAA - 7 aa", "# Database: C:\\Users\\Sebastian\\Downloads\\uniprot-taxonomy-lambda.fasta",
			"# 0 hits found", "# SSEARCH 36.3.8g Dec, 2017", "# Query: BBBBBB - 6 aa",
			"# Database: C:\\Users\\Sebastian\\Downloads\\uniprot-taxonomy-lambda.fasta",
			"# 0 hits found", "# SSEARCH 36.3.8g Dec, 2017", "# Query: IARVRDIKPVWALANDMNCSAG - 22 aa",
			"# Database: C:\\Users\\Sebastian\\Downloads\\uniprot-taxonomy-lambda.fasta",
			"# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score",
			"# 2 hits found", "IARVRDIKPVWALANDMNCSAG\tsp|P03711|SCAF_LAMBD\t100.00\t22\t0\t0\t1\t22\t147\t168\t1e-012\t59.5",
			"IARVRDIKPVWALANDMNCSAG\tsp|P36273|VG05_BPP21\t66.67\t21\t7\t0\t1\t21\t143\t163\t2.7e-006\t38.4",
			"# SSEARCH processed 3 queries")
		parse_ssearch_out(sample_out)

	}, "data.frame")

})

test_that("Generating temporary fasta files work", {
	seqs <- c("AAAAAAA", "BBBBBB", "IARVRDIKPVWALANDMNCSAG")
	foo <- make_tmp_fasta(seqs)

	expect_true(file.exists(foo))
	expect_equal(length(readLines(foo)), 6)
})


test_that("Wrapper arround ssearch works", {
	skip_if(system('ssearch') == "127", 'ssearch not in path, skipping wrapper testing')
	seqs <- c("AAAAAAA", "BBBBBB", "IARVRDIKPVWALANDMNCSAG")
	foo <- find_in_fasta(
		sequences = seqs,
		target_fasta = 	"C:\\Users\\Sebastian\\Downloads\\uniprot-taxonomy-lambda.fasta",
		ssearch_exe = "C:\\Users\\Sebastian\\opt\\fasta36\\ssearch36.exe")
})



