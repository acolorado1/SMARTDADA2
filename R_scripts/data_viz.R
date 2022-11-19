# needed libraries 
require(ggplot2)
require(reshape2)
require(dplyr)
 

#' Boxplot of quality scores per positions 
#'
#' @param dataframe A dataframe where colnames are positions on the read and 
#' rows contain the quality score found at that position per read.
#' @param threshold Integer default value 20. Boxes will be colored by whether 
#' their average quality score is above or below the threshold. 
#'
#' @return Boxplot as png in working directory 
#' @export
#'
#' @examples distribution_boxplot(df, 20)
#' 
distribution_boxplox <- function(dataframe, 
                                 threshold = 20){
  
  meltdata <- melt(dataframe)
  meltdata$variable <- as.factor(meltdata$variable)
  meltdata <- meltdata %>%
    group_by(variable)%>% 
    mutate(Average = ifelse(mean(value) >= theshold, 
                            as.character(threshold), 
                            as.character(threshold)))
  
  ggplot(meltdata, aes(x = variable, y = value, fill = Average)) + 
    geom_boxplot()+ 
    xlab("Position") + 
    ylab("Quality Scores")+ 
    ggtitle("Distribution of Quality Scores by \n Position on Read")+
    theme_bw()+ 
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave("plots/boxplot.png")
}




#' Line plot of average expected error per position 
#'
#' @param df Dataframe containing two columns: Position and average expected 
#' error at the position. 
#'
#' @return Save line plot to working directory 
#' @export
#'
#' @examples average_error_lineplot(df)
#' 
average_error_lineplot <- function(df){ 
  df$Position <- as.factor(df$Position)
  ggplot(df, aes(x = Position, y = AverageExpectedError, group = 1))+ 
    geom_line()+
    ylab("Average Expected Error")+ 
    ggtitle("Average Expected Error by Position")+
    theme_bw()+ 
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave('plots/LinePlotAvgEE.png')
}


