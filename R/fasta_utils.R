
#' Make a temporary fasta file
#'
#' @param sequences Character vector of sequences to be written
#' @param seq_names Character vector of names to be used in the fasta header (defaults to the sequences)
#'
#' @return Path to the temporary fasta
#' @export
#'
#' @examples
#' make_tmp_fasta(c("PEPTIDE", "TIDALES"))
#' make_tmp_fasta(c("PEPTIDE", "TIDALES"), seq_names = c("Pep1", "Pep2"))
make_tmp_fasta <- function(sequences, seq_names = NULL ) {
	fn  <- tempfile(pattern = "", fileext = ".fasta")
	print(paste0("Writting fasta to ", fn))
	if (is.null(seq_names)) {
		seq_names <- sequences
	}

	zz <- file(fn, "w")
	for (i in seq_along(sequences)) {
		writeLines(paste0(">", seq_names[i]), con = zz)
		writeLines(sequences[i], con = zz)
	}
	close(zz)
	return(fn)
}




#' R wrapper to ssearch
#'
#' Generates a temporary fasta file with the sequences provided and returns the
#' matches in the provided fasta file.
#'
#' It wraps sseach to do the search and then returns the console output
#'
#' @param sequences sequences to be searched
#' @param target_fasta path to fasta file containing the target sequences
#' @param ssearch_exe location of the ssearch executable (or alias in your path)
#' @param seq_names optional names to be given to the sequences
#' @param min_eval minimum E value to be reported, defaults to 0.01
#'
#' @return tabular console output of ssearch as a character vector (element per line)
#' @export
#'
#' @examples
#' \donttest{find_in_fasta(c('PEPTIDE'), 'arabidopsis.fasta', '~/usr/bin/ssearch36')}
find_in_fasta <- function(
	sequences,
	target_fasta,
	ssearch_exe,
	seq_names = NULL,
	min_eval = 0.01) {
	# Returns a list of DFs
	tmpfasta <- make_tmp_fasta(sequences = sequences, seq_names = seq_names)
	out <- system2(ssearch_exe,
				   args = paste("-m 8C -E",min_eval , tmpfasta, target_fasta),
				   stdout = TRUE)
	return(out)

}


#' Parses the console output of ssearch
#'
#' It is designed to parse que tabular out of ssearch (flag -m 8C)
#'
#' @param ssearch_out tabular console output of ssearch as a character vector (element per line)
#'
#' @return data.frame with output of the ssearch search
#' @export
#'
#' @examples
#' \donttest{
#'     ss_out <- find_in_fasta(c('PEPTIDE'), 'arabidopsis.fasta', '~/usr/bin/ssearch36')
#'     parse_ssearch_out(ss_out)
#'     }

parse_ssearch_out <- function(ssearch_out) {
	# Output from "-m 8C" flag

	RE <- ".* Fields: (.*)($|\\n)"
	col_names_line <- ssearch_out[which(grepl(RE, ssearch_out))[[1]]]
	col_names <- unlist(strsplit(gsub(RE, "\\1", col_names_line), ','))
	if (length(ssearch_out) > 1) {
		ssearch_out <- paste0(ssearch_out, collapse = "\n")
	}
	data.table::fread(
		ssearch_out,
		header = FALSE,
		col.names = col_names)

}
