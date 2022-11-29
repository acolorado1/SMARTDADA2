require(readr)
require(dplyr)
require(reshape2)
require(ggplot2)
require(argparse)
require(gt)
require(webshot2)

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

ggsave("plots/ReadLengthByAvgEE.png",
       width = 10,
       height = 10)

dev.off()

# line plot of read counts over and under the maxEE
melt_param_info <- melt(parameter_info, id.vars = 'ReadLength')
melt_param_info <- melt_param_info %>% 
  filter(variable == "ReadsUnderMaxEE" | variable == "ReadsOverMaxEE")

ggplot(melt_param_info, aes(x = ReadLength, y = value)) + 
  geom_line(aes(color = variable)) +
  scale_colour_manual(values=c(ReadsUnderMaxEE="#00FF00",ReadsOverMaxEE="#FF0000"))+
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  xlab('Read Length (bp)') + 
  ylab('Number of reads') + 
  guides(color=guide_legend(title=" "))

ggsave("plots/ReadCountOverMaxEE.png",
       width = 10,
       height = 10)

dev.off()


# creating a table showing read length and 25% quantile of AvgEEPerPosition
quantiles <- quantile(parameter_info$AvgEEPerPosition)

summary_table <- parameter_info %>% 
  filter(ReadLength == max(ReadLength) | AvgEEPerPosition <= quantiles[[2]]) %>%
  arrange(-ReadLength, AvgEEPerPosition) %>% 
  select(Indexes, ReadLength, AvgEEPerPosition) %>%
  distinct()


max_readlen <-  max(parameter_info$ReadLength)
min_avgEE <- min(parameter_info$AvgEEPerPosition)

gt_summary <- summary_table %>% 
  gt() %>% 
  tab_style(
    style = cell_fill(color = "gray65"),
    locations = cells_body(rows = ReadLength == max_readlen)
  )%>% 
  tab_style(
    style = cell_fill(color = "gray65"), 
    locations = cells_body(rows = AvgEEPerPosition == min_avgEE)
    ) %>% 
  gtsave(filename = "plots/t_len_by_AvgEE.png")


# creating a table of Read length and 75% quantile of ReadsOverMaxEE
quantiles <- quantile(parameter_info$ReadsOverMaxEE)

summary_table <- parameter_info %>% 
  filter(ReadsOverMaxEE <= quantiles[[4]]) %>% 
  arrange(-ReadLength, AvgEEPerPosition) %>% 
  distinct()

gt_summary <- summary_table %>% 
  gt() %>% 
  tab_style(
    style = cell_fill(color = "gray65"),
    locations = cells_body(rows = ReadLength == max_readlen)
  )%>% 
  tab_style(
    style = cell_fill(color = "gray65"), 
    locations = cells_body(rows = AvgEEPerPosition == min_avgEE)
  ) %>% 
  gtsave(filename = "plots/t_len_by_ReadsOverMaxEE.png")



