
#' Make a temporary fasta file
#'
#' @param sequences Character vector of sequences to be written
#' @param seq_names Character vector of names to be used in the fasta header (defaults to the sequences)
#' @param warn_unique_names bool, logical indicating if the function should check that any names are duplicated
#'
#' @return Path to the temporary fasta
#' @export
#'
#' @examples
#' make_tmp_fasta(c("PEPTIDE", "TIDALES"))
#' make_tmp_fasta(c("PEPTIDE", "TIDALES"), seq_names = c("Pep1", "Pep2"))
make_tmp_fasta <- function(sequences, seq_names = NULL, warn_unique_names = TRUE ) {
	fn  <- tempfile(pattern = "", fileext = ".fasta")
	message(paste0("Writting fasta to ", fn))
	if (is.null(seq_names)) {
		seq_names <- sequences
	}

	if (warn_unique_names) {
		if (any(duplicated(sequences))) { warning('Not all names are unique') }
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
#' @param add_flags additional flags to be passed to ssearch
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
	min_eval = 0.01,
	add_flags = '') {
	# Returns a list of DFs
	tmpfasta <- make_tmp_fasta(
		sequences = sequences,
		seq_names = seq_names)
	out <- system2(
		ssearch_exe,
		args = paste(
			"-m 8C -E",
			min_eval,
			add_flags,
			tmpfasta,
			target_fasta),
		stdout = TRUE)
	return(out)

}


#' Find Exact matches for sequences in a fasta file
#'
#' Its a simple wrapper arround `Biostrings::readAAStringSet` and
#'  `Biostrings::vcountPattern`.
#'
#' @param sequences string vector
#' @param target_fasta file location
#' @param in_mem_fasta output of Biostrings::readAAStringSet stored in memmory
#'
#' @return number of fasta entries in target that match each sequence
#' @export
#'
#' @examples
find_exact_in_fasta <- function(
	sequences,
	target_fasta,
	in_mem_fasta = NA) {
	if (!requireNamespace("Biostrings", quietly = TRUE)) {
		stop("Biostrings (Bioconductor) is not installed, please install to use this function")
	}

	if(!is.na(in_mem_fasta)) {
		targets <- in_mem_fasta
	} else {
		targets <- Biostrings::readAAStringSet(target_fasta)
	}
	out <- purrr::map_int(sequences, .f = (
		function(x) {
			sum(as.logical(Biostrings::vcountPattern(x, targets)))
		}
	))
	return(out)
}


#' Parses the console output of ssearch
#'
#' It is designed to parse que tabular output of ssearch (flag -m 8C)
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
#' \donttest{
#'     parse_ssearch_out(readLines('./ssearch.out.txt'))
#'     }
parse_ssearch_out <- function(ssearch_out) {
	# Output from "-m 8C" flag

	RE <- ".* Fields: (.*)($|\\n)"
	col_names_line <- ssearch_out[which(grepl(RE, ssearch_out))[[1]]]
	col_names <- unlist(strsplit(gsub(RE, "\\1", col_names_line), ','))
	ssearch_out <- ssearch_out[!grepl('^#', ssearch_out)]
	if (length(ssearch_out) > 1) {
		ssearch_out <- paste0(ssearch_out, collapse = "\n")
	}
	data.table::fread(
		ssearch_out,
		header = FALSE,
		col.names = col_names)

}



#' Query the current files in the uniprot FTP
#'
#' Queries the uniprot ftp server and returns a data frame of the
#' served files for a given domain (Eukaryota, Bacteria, Archaea or Viruses)
#'
#' @param domain string. any of Eukaryota, Bacteria, Archaea or Viruses
#' @param dry logical. If false no files will be downloaded and will print what would have.
#'
#' @return data frame with the files served under that domain
#' @export
#'
#' @examples
query_current_uniprot <- function(domain = 'Eukaryota', dry = FALSE) {
	tmpfile <- tempfile(fileext = '.txt')
	url_query <- glue::glue(
		'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/',
		'knowledgebase/reference_proteomes/',
		'{domain}/',
		domain = domain)
	if (dry) {
		warning('Dry run, not downloading anything, printing URLS')
		print(url_query)
		return(FALSE)
	} else {
		print("Querying Uniprot: ...")
		curl::curl_download(
			url_query,
			tmpfile)
		return(data.table::fread(tmpfile))
	}


}


#' Downloads from uniprot the current release of given proteomes
#'
#' Downlaods using the ftp client from unioprot the given fasta files that
#' correspond to the proteome or species id provided
#'
#' @param proteome_id For example 64320 for the zika virus
#' @param species_id For example 9606 for human
#' @param domain any of Eukaryota, Bacteria, Archaea or Viruses
#' @param prefix string. A prefix to be given to the downloaded files
#' @param dest_dir string. Valid path of a directory where the file will be downloaded, will be created if inexistent
#' @param dry logical. If false no files will be downloaded and will print what would have.
#' @param force logical. If FALSE will check if the files exist and will not overwrite them.
#' @param cached_uniprot_release data.frame. the output of `query_current_uniprot`, when existent it will skip that query and use the provided one (accelerates drastically the download of multiple proteomes)
#' @param ungz logical. If TRUE the files downloaded as .gz will be uncompressed in the directory after downloading them
#'
#' @return TRUE when completed succesfully, false when it is a dry run
#' @export
#'
#' @examples
download_ref_proteome_fasta <- function(
	proteome_id = NULL,
	species_id = NULL,
	domain = 'Eukaryota',
	prefix = '',
	dest_dir = '.',
	dry = TRUE,
	force = TRUE,
	cached_uniprot_release = NULL,
	ungz = TRUE) {

	filter_args <- c(proteome_id, species_id)

	if (all(is.null(filter_args))) {
		stop('Please Provide a valid proteome or a species id')
	}

	# Check destination directory as valid
	if (!dir.exists(dest_dir)) {
		tryCatch({
			dir.create(dest_dir)
		},
		error = function(e){
			print(e)
			stop(glue::glue(
				'Destination directory does not exist and could not create it,',
				' check that a valid path is being provided'))
		})
	}

	if (!is.null(cached_uniprot_release)) {
		current_UP <- cached_uniprot_release

	} else if (!dry) {
		current_UP <- query_current_uniprot(domain = domain)
	}

	# TODO make a subfunction from here to be able to pass a list of proteomes
	# or proteome ids

	REGEX = glue::glue(
		"{proteome_id}_{species_id}(_additional)?.fasta.gz",
		proteome_id = ifelse(is.null(proteome_id), '.*', proteome_id),
		species_id = ifelse(is.null(species_id), '.*', species_id))

	# V9 is the filed where files names lay in the current uniprtot release
	downloadable_files <- grep(
		REGEX,
		current_UP[['V9']],
		value = TRUE)


	filenames <- Vectorize(function(filename) {
		glue::glue(
			'{dest_dir}/{prefix}{downloadable_file}',
			dest_dir = dest_dir,
			prefix = prefix,
			downloadable_file = filename)
	}, 'filename', USE.NAMES = FALSE) (downloadable_files)

	urls <- Vectorize(function(target) { glue::glue(
		'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/',
		'knowledgebase/reference_proteomes/',
		'{domain}/{target}',
		domain = domain,
		target = target)
	}, 'target', USE.NAMES = FALSE) (downloadable_files)

	if (dry) {
		warning('Dry run, not downloading anything, printing URLS')
		cat("\n urls to be downloaded \n ===============\n")
		print(urls)
		cat("\n File names to be used \n ===============\n")
		print(filenames)
		return(FALSE)
	} else {
		for (i in seq_along(urls)){
			print(
				glue::glue(
					"Downloading {url}, to {filename}",
					url = urls[i],
					filename = filenames[i]))
			file.create(filenames[i])
			curl::curl_download(url = urls[i], destfile = filenames[i])
			if (ungz) {
				stopifnot(grepl("(.*\\.\\w+).gz", filenames[i]))
				writeLines(
					readLines(
						filenames[i]),
					gsub(
						"(.*\\.\\w+).gz",
						"\\1",
						filenames[i]))
			}
			return(TRUE)
		}
	}
}




#' Get mass difference between two peptides
#'
#' @param base str, base peptide sequence to be compared (diff will be positive if this one is bigger).
#' @param replacement str, peptide sequence to be compared against.
#' @param verbose bool, indicates wether print and message statements should be used, false indicades it should be silent.
#'
#' @return double, mass difference betweent he two peptides
#' @export
#'
#' @examples
get_massdiff <- function(base, replacement, verbose = FALSE) {
	if (any(is.na(c(base, replacement))) & verbose) {
		if (verbose) warning('Either Base or replacement sequenes are NA')
	}
	suppressWarnings({
		massdiff <- Peptides::mw(replacement, monoisotopic = TRUE) -
			Peptides::mw(base, monoisotopic = TRUE)
	})

	if (verbose) {
		print(
			glue::glue(
				'{first} --> {second} = {massdiff}',
				first = base,
				second = replacement,
				massdiff = massdiff)
		)}
	return(massdiff)
}



#' Calculates mass differenes by taking a cigar string as an input
#'
#' @param str string, cigar string (or vector of cigar strings). ie: '4QK3SE1TY1KR'
#' @param verbose bool, indicates wether print and message statements should be used, false indicades it should be silent.
#'
#' @return double, vector of mass differences
#' @export
#'
#' @examples
cigar_massdiff <- Vectorize(function(str, verbose = FALSE) {
	replacement_pairs <- stringi::stri_extract_all(str, regex = '[A-Z]{2}')
	replacement_pairs <- plyr::laply(replacement_pairs, strsplit, '')
	foo <- lapply(
		replacement_pairs,
		function(x){
			tmp <- get_massdiff(
				x[1], x[2],
				verbose = verbose)
			return(sum(tmp))
		})
	return(sum(unlist(foo)))
}, vectorize.args = 'str', USE.NAMES = FALSE)







