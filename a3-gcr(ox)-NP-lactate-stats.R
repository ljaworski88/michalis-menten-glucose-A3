# This is a script used to analyze any differences between Lactate Ratios in NP
# cells cultured under differing oxygen tensions over a 24 hour adaptation
# period

# Author <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2017
###############################################################################
## Library imports
###############################################################################
library(ggplot2)
###############################################################################
## Functions
###############################################################################
## Summarizes data.This function was copied from cookbook-r.com
## Gives count, mean, standard deviation, standard error of the mean, and 
## confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Used for comparison, add lactate production rate -> get glucose consumption rate
FromLactate <- function(lactate.production.rate){
  glucose.consumption.rate <- lactate.production.rate/
    lactate.summary.avg[['lactate.ratio']]
  return(glucose.consumption.rate)
}

###############################################################################
## Script
###############################################################################
# Retrieving and arranging the data

# The directory where the datafile is kept
if (Sys.info()[['sysname']]=='Linux'){
  working.directory <- paste0('/home/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/Lactate_Assay/',
                              'a3-gcr(ox)-NP/Data/')
} else {
  working.directory <- paste0('C:/Users/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/Lactate_Assay/',
                              'a3-gcr(ox)-NP/Data/')
}

setwd(working.directory)

# the lactate ratio numbers, arranged in a csv file by copying the lactate ratio numbers
# from the 'Summary p19-p23.xlsx' file, in the sheet final numbers
raw.numbers <- read.csv('lactate-summary-p19_p23.csv', header = FALSE)

# assigning names to the data columns: sName-Sample Name; lactate.ratio-Lactate Ratio
names(raw.numbers)[1] <- 'sName'
names(raw.numbers)[2] <- 'lactate.ratio'

# changing the sample names to characters instead of factors
raw.numbers$sName<-as.character(raw.numbers$sName)

# adding a column to identify from which pig each sample came from
# this will be the 'blocking' variable in the ANOVA
raw.numbers$pig <- NA
raw.numbers$pig[grep('P19', raw.numbers$sName)] <- 'P19'
raw.numbers$pig[grep('P20', raw.numbers$sName)] <- 'P20'
raw.numbers$pig[grep('P21', raw.numbers$sName)] <- 'P21'
raw.numbers$pig[grep('P22', raw.numbers$sName)] <- 'P22'
raw.numbers$pig[grep('P23', raw.numbers$sName)] <- 'P23'
raw.numbers$pig <- as.factor(raw.numbers$pig)

# adding a column identifying the experimental groups
# this is our independant variable
raw.numbers$O2group <- NA
raw.numbers$O2group[grep('-2.5%', raw.numbers$sName)] <- '2.5%'
raw.numbers$O2group[grep('-5%', raw.numbers$sName)]   <- '5%'
raw.numbers$O2group[grep('-21%', raw.numbers$sName)]  <- '21%'
raw.numbers$O2group <- factor(raw.numbers$O2group)

# arrangeing the order of the experimental groups so that when they
# are graphed later they are displayed in the correct order
levels(raw.numbers$O2group)<-c('2.5%', '5%', '21%')

# running the ANOVA and the post hoc tests as well as displaying the results.
lactate.anova <- aov(lactate.ratio ~ O2group + pig, data = raw.numbers)
summary.aov(lactate.anova)
TukeyHSD(lactate.anova)

# getting the descriptive stats so that the standard error
# can be used for the error bars
lactate.summary     <- summarySE(raw.numbers, 
                              measurevar="lactate.ratio",
                              groupvars="O2group")
lactate.summary.avg <- summarySE(raw.numbers, 
                              measurevar="lactate.ratio")

###############################################################################
# Graphng and displaying the results

# making a bar graph of the lactate ratios
lactate.graph <- ggplot(lactate.summary, 
                        aes(O2group, lactate.ratio, fill=O2group)) +
  xlab('Oxygen Tension') + 
  ylab('LPR:GCR') +
  geom_errorbar(aes(ymin=lactate.ratio,
                    ymax=lactate.ratio+se),
                width = .2) +
  geom_col(color='black') +
  # ggtitle('Lactate Production Ratio Vs\nCultured Oxygen Tension') +
  scale_fill_manual(values=c('white', 'grey', 'black')) +
  theme_bw() +
  theme(legend.position='none')

plot(lactate.graph)