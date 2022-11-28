# install.R install all required packages

# read a csv file that contains the package name and the version
# -- make the first column as the sequence of required packages
# -- -- if not installed, use

RDepInstall <- function() {

    requiredLibs <- list("BiocManager", "ggplot2", "reshape2", "dplyr")
    for (lib in requiredLibs) {

        # checking if package are installed if not they are installed
        if (!require(lib, quietly = TRUE))
            # print("Installing: %s", lib)
            install.packages(lib, repos='http://cran.us.r-project.org', quiet=TRUE)

    }
}