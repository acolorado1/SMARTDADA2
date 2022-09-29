# needed libraries 
library(ggplot2)
library(reshape2)
library(dplyr)

# plot distributions per position 
distribution_boxblox <- function(dataframe){
  
  meltdata <- melt(dataframe)
  meltdata$variable <- as.factor(meltdata$variable)
  meltdata <- meltdata %>%
    group_by(variable)%>% 
    mutate(Average = ifelse(mean(value) >= 20, '>20', '<20'))
  
  ggplot(meltdata, aes(x = variable, y = value, fill = Average)) + 
    geom_boxplot()+ 
    xlab("Position") + 
    ylab("Quality Scores")+ 
    ggtitle("Distribution of Quality Scores by \n Position on Read")+
    theme_bw()+ 
    theme(plot.title = element_text(hjust = 0.5))
}




# plot line plot of expected errors 
average_error_lineplot <- function(df){ 
  df$Position <- as.factor(df$Position)
  ggplot(df, aes(x = Position, y = AverageExpectedError, group = 1))+ 
    geom_line()+
    ylab("Average Expected Error")+ 
    ggtitle("Average Expected Error by Position")+
    theme_bw()+ 
    theme(plot.title = element_text(hjust = 0.5))
}


