require(readr)
require(dplyr)
require(reshape2)
require(ggplot2)
require(argparse)

parser <- ArgumentParser()
parser$add_argument('--tsv_file', help = "tsv file containing parameter info")
args <- parser$parse_args()


parameter_info <- read_delim(args$tsv_file,
                             delim = "\t", escape_double = FALSE,
                             col_types = cols(...1 = col_integer(),
                                              Indexes = col_character(),
                                              ReadLength = col_integer(),
                                              ReadsUnderMaxEE = col_integer(),
                                              ReadsOverMaxEE = col_integer()),
                             trim_ws = TRUE)


# scatter plot of average expected error per position by read length 
ggplot(parameter_info, aes(x = ReadLength, y = AvgEEPerPosition)) + 
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ylab('Average Expected Error Per Position') + 
  xlab('Read Length (bp)') 

ggsave("ReadLengthByAvgEE.png")


# line plot of read counts over and under the maxEE
melt_param_info <- melt(parameter_info, id.vars = 'ReadLength')
melt_param_info <- melt_param_info %>% 
  filter(variable == "ReadsUnderMaxEE" | variable == "ReadsOverMaxEE")

ggplot(melt_param_info, aes(x = ReadLength, y = value)) + 
  geom_line(aes(color = variable)) +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  xlab('Read Length (bp)') + 
  ylab('Number of reads') + 
  guides(color=guide_legend(title=" "))

ggsave("ReadCountOverMaxEE.png")
