
# This series of scripts rely on proteowizzard being installed
# and in your path.


parse_msaccess_out <- function(
	file,
	force_R_filenames = FALSE,
	append_metadata = TRUE) {

	# col_names <- scan(file, what = character(), nlines = 1, skip = 1)
	DT <- data.table::fread(file)

	if (force_R_filenames) {
		data.table::setnames(
			DT,
			colnames(DT),
			make.names(colnames(DT), unique = TRUE))
	}

	if (append_metadata) {
		# Removes path from the file provided ...
		# might have a bug when escape characters are provided
		filedata <- gsub(".*[/\\](.*).tsv", '\\1', file)
		filedata <- regmatches(
			filedata,
			regexec(
				'(.*).(?:tic|sic).([\\d+-.]+)',
				filedata,
				perl = TRUE))[[1]]
		names(filedata) <- c('datafile', 'raw_file', 'mass_range')

		for (i in seq_along(filedata)) {
			DT[[names(filedata)[i]]] <- filedata[i]
		}
	}
	return(DT)
}



# Get XIC
get_xic <- function(
	raw_files,
	mz = NA,
	mz_width = NA,
	mz_max = mz + (mz_width/2),
	mz_min = mz - (mz_width/2),
	num_isotopes = 2,
	charge = 2,
	rt = NA,
	rt_window = NA,
	rt_start = rt - (rt_window/2),
	rt_end = rt + (rt_window/2),
	out_folder = NA,
	msaccess_executable = 'msaccess',
	ms_level = 1,
	dry = TRUE) {

	if (!dry) {
		if (is.na(out_folder)) {
			out_folder <- tempdir()
			msg <- glue::glue(
				"Created directory in {dir} to store the output of the command",
				dir = out_folder)
			message(msg)
		}

		if (!dir.exists(out_folder)) {
			out_folder <- file.path()
			dir.create(out_folder)
		}
	}

    numeric_args <- c(unlist(mz_max), unlist(mz_min), rt_start, rt_end)

	if (!all(is.numeric(numeric_args)) |
		any(is.na(numeric_args))) {
		errormsg <- glue::glue(
			'Arguments needed for the command are not numeric, check input',
			'mz_min = {mz_min},',
			'mz_max = {mz_max},',
			'rt_start = {rt_start},',
			'rt_end = {rt_end},'
		)
		stop(errormsg)
		# TODO Find a better way to assert the validity of arguments
	}

    isotope_mods <- sort(
    	rep(seq_len(num_isotopes)-1) * (1/charge), length(mz_min))


	execution_commands <- glue::glue(
		'-x "tic mz=[{mz_min},{mz_max}] delimiter=tab"',
		mz_min = mz_min + isotope_mods,
		mz_max = mz_max + isotope_mods
	)
	print(execution_commands)

	call <- glue::glue(
		'{raw_files} {execution_command}',
		' --outdir {out_folder}',
		' --filter="msLevel {ms_level}" ',
		' --filter "scanTime [{rt_start}, {rt_end}]" -v',
		raw_files = paste(
			glue::glue('"{raw_files}"',
				 raw_files = raw_files),
			collapse = ' '),
		execution_command = paste(execution_commands, collapse = ' '),
		ms_level = ms_level,
		rt_start = rt_start,
		rt_end = rt_end
	)

	# Assert that msaccess is in the path
	ifelse(
		dry,
		return(paste(msaccess_executable, call)),
		stout <- system2(msaccess_executable, call, stdout = TRUE)
		)
	message('Output from the external command')
	message(paste(stout), collapse = '\n')

	out_files <- stout[grepl('Writing file', stout)]
	out_files <- gsub("^.*Writing file (.*)$", '\\1', out_files)

	out_DF <- purrr::map_dfr(out_files, parse_msaccess_out)

	return(out_DF)
}

# This would be how to plot the output
# foo %>% ggplot2::ggplot(ggplot2::aes(x =  rt, y = sumIntensity, colour = mass_range)) + ggplot2::geom_line() + ggplot2::scale_y_sqrt() + theme_bw()


# Get SIM
#
#
	call <- glue::glue(
'msaccess.exe {raw_files} -x "sic mzCenter={mz} radius={radius}',
'radiusUnits=amu delimiter=tab" ',
'--filter="msLevel {1}" --filter "scanTime [{rt_start},{rt_end}]" -v'
)
