# Check dependecies 
packages <- c("readr", "dplyr", "reshape2", "ggplot2", "argparse")

for (package in packages) {
  if(!require(package, character.only = T)){
    install.packages(package)
    library(package)
  }
}

# create parser 
parser <- ArgumentParser()
parser$add_argument('--trim_file', help = "TSV file containing trim info")
parser$add_argument('--sumEE_file', help = 'TSV file containing sum of EE info')
args <- parser$parse_args()


# Read in data 
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
## Scatter plot
sp <- ggplot(trim_info, aes(ReadLength, AvgEEPerPosition)) +
  geom_point() +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ylab('Average Expected Error Per Position') +
  xlab('Read Length (bp)')

ggsave("output/plots/ScatterReadLengthByAvgEE.png",
       plot = sp,
       width = 10,
       height = 10)

## Heatmap
lefttrim <- unique(trim_info$LeftTrim)
righttrim <- unique(trim_info$RightTrim)
hm <- ggplot(trim_info, aes(LeftIndex, RightIndex, fill= AvgEEPerPosition)) +
  geom_tile() +
  geom_vline(xintercept = lefttrim) + 
  geom_hline(yintercept = righttrim)+
  scale_y_continuous(trans = "reverse", breaks = unique(trim_info$RightIndex)) +
  scale_fill_gradient(low="white", high="blue") +
  theme_bw()

ggsave("output/plots/HeatmapIndexValueByAvgEE.png",
       plot = hm,
       width = 10,
       height = 10)




# read quality by retained reads
## Histogram
melt_sumEEInfo <- melt(sumEE_info)

hi <- ggplot(melt_sumEEInfo, aes(x = value, fill = variable)) +
  geom_histogram(position = "dodge") +
  geom_vline(xintercept = 2) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(limits = c(0,4)) +
  guides(fill=guide_legend(title=" "))+
  xlab('Sum of Expected Error')


ggsave("output/plots/HistogramRetainedReadCount.png",
       plot = hi,
      width = 10,
      height = 10)



