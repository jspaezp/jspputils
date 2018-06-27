context("uniprot_querying")

tmpdir <- tempdir()
UP_db_virus <- NULL
UP_db_euk <- NULL
UP_db_bacteria <- NULL

test_that("Can Download the list of served reference proteomes", {
	expect_output({
		str(query_current_uniprot())
	}, 'data.frame')
	expect_output({
		UP_db_virus <<- query_current_uniprot(domain = 'Viruses')
		str(UP_db_virus)
	}, 'data.frame')
	expect_output({
		UP_db_euk <<- query_current_uniprot(domain = 'Eukaryota')
		str(UP_db_euk)
	}, 'data.frame')
	expect_output({
		UP_db_bacteria <<- query_current_uniprot(domain = 'Bacteria')
		str(UP_db_bacteria)
	}, 'data.frame')
})


test_that("An error is thrown when invalid arguments are given for proteome download",{
	expect_error(download_ref_proteome_fasta())
})

test_that("Dry run returns false and does not call curl", {
	expect_false({
		download_ref_proteome_fasta(
			proteome_id = 'UP000054557',
			species_id = '64320',
			prefix = 'Zika',
			domain = 'Viruses',
			dry = TRUE,
			dest_dir = tmpdir,
			cached_uniprot_release = UP_db_virus)
	})
})

skip_if(
	is.null(UP_db_virus),
	'Was not able to get the uniprot DB, skipping depending tests'
	)


test_that("Can Download proteomes by proteomeID", {
	expect_true({
		download_ref_proteome_fasta(
			proteome_id = 'UP000054557',
			species_id = '64320',
			prefix = 'Zika',
			domain = 'Viruses',
			dry = FALSE,
			dest_dir = tmpdir,
			cached_uniprot_release = UP_db_virus)
	})
})


test_that("Can Download proteomes by organismID", {
	expect_true({
		download_ref_proteome_fasta(
			species_id = '308745',
			prefix = 'Aspergillus_rambellii',
			dry = FALSE,
			dest_dir = tmpdir,
			cached_uniprot_release = UP_db_euk)
	})
	expect_true({
		download_ref_proteome_fasta(
			species_id = '64320',
			prefix = 'Viruses',
			domain = 'Viruses',
			dry = FALSE,
			dest_dir = tmpdir,
			cached_uniprot_release = UP_db_virus)
	})
})



test_that("Can Download proteomes by proteomeID and organismID", {
	expect_true({
		download_ref_proteome_fasta(
			proteome_id = 'UP000054557',
			species_id = '64320',
			prefix = 'Viruses',
			domain = 'Viruses',
			dry = FALSE,
			dest_dir = tmpdir,
			cached_uniprot_release = UP_db_virus)
	})
})
