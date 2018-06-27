
# This series of scripts rely on proteowizzard being installed
# and in your path.


# Get XIC
get_xic <- function(
	raw_file,
	mz = NULL,
	mz_width = NULL,
	mz_max = mz + (mz_width/2),
	mz_min = mz - (mz_width/2),
	num_isotopes = 2,
	rt = NULL,
	rt_window = NULL,
	rt_start = rt - (rt_window/2),
	rt_end = rt + (rt_window/2),
	out_folder = NULL) {

	print(c(mz_max, mz_min, rt_start, rt_end))

	if (is.null(out_folder)) { out_folder <- tempdir() }

	if (!dir.exists(out_folder)) {
		out_folder <- file.path()
		dir.create(out_folder)
	}

	if (all(is.numeric(...))) { # Find a better way to assert the validity of arguments
		# Assert that msaccess is in the path
		call <- glue::glue('')
		system2(call, )
	}

	return(TRUE)
}

# Get SIM
