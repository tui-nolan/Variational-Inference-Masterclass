############### R-function: construct_txt ###############

# For combining a set of files into a joint txt file

# Created: 16 AUG 2022
# Last changed: 16 AUG 2022

construct_txt <- function(directory = ".", pattern = ".txt", file_name, header) {
	
	files <- list.files(pattern = pattern)
	file_numbers <- as.numeric(regmatches(files, regexpr("[0-9]+", files)))
	files <- files[order(file_numbers)]
	n_files <- length(files)
	
	col_names <- names(read.table(files[[1]], header = header))
	n_col <- length(col_names)
	write(col_names, file_name, ncol = n_col, append = FALSE)
	
	for(i in 1:n_files) {
		
		res <- unname(unlist(read.table(files[[i]], header = header)))
		write(res, file_name, ncol = n_col, append = TRUE)
	}
}

############### End construct_txt ###############