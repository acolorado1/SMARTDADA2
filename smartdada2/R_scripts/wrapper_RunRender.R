packages <- c("argparse")

for (package in packages) {
  if(!require(package, character.only = T)){
    install.packages(package)
    library(package)
  }
}

parser <- ArgumentParser()
parser$add_argument('--trim_file', help = "TSV file containing trim info")
parser$add_argument('--sumEE_file', help = 'TSV file containing sum of EE info')
args <- parser$parse_args()

render_report <- function(trim_file, sumEE_file) {
  rmarkdown::render(
    "smartdada2/R_scripts/SMARTDADA2_InteractiveOutput.Rmd", params = list(
      trim_file = trim_file,
      sumEE_file = sumEE_file
    ),
    output_file = paste0("../../output/SMARTDADA2_InteractiveOutput.html")
  )
}

render_report(args$trim_file, args$sumEE_file)