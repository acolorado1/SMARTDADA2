require(readr)
require(dplyr)
require(reshape2)
require(ggplot2)
require(argparse)

parser <- ArgumentParser()
parser$add_argument('--trim_file', help = "TSV file containing trim info")
parser$add_argument('--sumEE_file', help = 'TSV file containing sum of EE info')
args <- parser$parse_args()


trim_info <- read_delim(args$trim_file, 
                        delim = "\t", escape_double = FALSE, 
                        col_types = cols(LeftIndex = col_integer(), 
                                         RightIndex = col_integer(), ReadLength = col_integer(), 
                                         ReadsUnderMaxEE = col_integer(), 
                                         ReadsOverMaxEE = col_integer()), 
                        trim_ws = TRUE)

sumEE_info <- read_delim(args$sumEE_file, 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

# plot read length vs quality 
ggplot(trim_info, aes(ReadLength, AvgEEPerPosition)) + 
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ylab('Average Expected Error Per Position') + 
  xlab('Read Length (bp)')

ggplot(trim_info, aes(LeftIndex, RightIndex, fill= AvgEEPerPosition)) + 
  geom_tile() + 
  scale_y_continuous(trans = "reverse", breaks = unique(parameter_info_LOZ_nano$RightIndex)) +
  scale_fill_gradient(low="white", high="blue") +
  theme_bw()

# read quality by retained reads
melt_sumEEInfo <- melt(sumEE_info)

ggplot(melt_sumEEInfo, aes(x = value, fill = variable)) + 
  geom_histogram(position = "dodge") + 
  theme_bw() + 
  theme(text = element_text(size = 20)) +
  scale_x_continuous(limits = c(0,4)) + 
  guides(fill=guide_legend(title=" "))
