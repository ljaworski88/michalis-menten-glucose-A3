# This is a script used to analyze any differences between Michalis-Menten 
# Constants in NP cells cultured under differing oxygen tensions over a 24 hour
# adaptation period

# Author      <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI          <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2017
###############################################################################
## Library imports
###############################################################################
library(ggplot2)
library(gridExtra)
library(cowplot)
###############################################################################
## Functions
###############################################################################
## Summarizes data. This function was copied from cookbook-r.com
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

# Copying the curve_fit function from the python program
MichalisMenten <- function(glucose.concentration, 
                           initial.concentration, 
                           Km, 
                           Vmax){
  rho = 1. / .0002  # 1 million cells/ 0.0002 L (200uL)
  Vmax <- Vmax
  return(Km / (rho * Vmax) * log(initial.concentration / glucose.concentration) + 
           (initial.concentration - glucose.concentration) / (Vmax * rho))
}

# making a function for the stat_function graphing for a cleaner look 
# (instead of including it in the lambda function)
kinetics <- function(glucose.concentration, vmax.vals, km.vals, i) {
  consumption.rate <- (vmax.vals[i, 'Vmax'])*glucose.concentration/
    (km.vals[i, 'Km'] + glucose.concentration)
  return(consumption.rate)
}

# Adds on the units to the output
GetConsumptionRate <- function(glucose.concentration, vmax.vals, km.vals, i){
  consumption.rate <- kinetics(glucose.concentration, vmax.vals, km.vals, i)
  return(paste(consumption.rate, '(nmol/million cells/hr)'))
}
###############################################################################
## Script
###############################################################################
# Retrieving and arranging the data

# Set the working directory to where the data is kept, make it workstation agnostic
if (Sys.info()[['sysname']] == 'Linux'){
  working.directory <- paste0('/home/',Sys.info()[['login']],
                              '/Dropbox/Lab/Experiments/GCR/',
                              'a3-gcr(ox)-NP/calculations/curvefit')
} else if (Sys.info()[['nodename']] == 'LUKAS-PC'){
  working.directory <- paste0('F:/Dropbox/Lab/Experiments/GCR/',
                              'a3-gcr(ox)-NP/calculations/curvefit')
} else {
  working.directory <- paste0('C:/Users/',Sys.info()[['user']],
                              '/Dropbox/Lab/Experiments/GCR/',
                              'a3-gcr(ox)-NP/calculations/curvefit')
}
setwd(working.directory)

# Getting the different groups
pig.ids <- c('P19', 'P20', 'P21', 'P22', 'P23')
experimental.groups <- c('2.5%', '5%', '21%')
for(id in pig.ids){
  for (group in experimental.groups){
    temp.df <- read.xlsx(paste0(id, '-GCR.xlsx'), sheetName = paste0(group, '-NP'), startRow = 2)
    temp.df <- temp.df[1:length(na.omit(temp.df[1,])),]
    temp.df <- temp.df[,c(1,((ncol(temp.df)-1)/2+2):ncol(temp.df))]
  }
}
michalis.constants <- read.csv('Michalis_Constants_Python.csv') 
sample.curvefit    <- read.csv('P22-5%-1-graph_data.csv')

# Getting the Michalis-Constants for the sample curve-fit graph
curve.km   <- subset(subset(michalis.constants, Pig == 'P22'), 
                     Experimental_Groups == '5%-NP')$Km[1]
curve.vmax <- subset(subset(michalis.constants, Pig == 'P22'), 
                     Experimental_Groups == '5%-NP')$Vmax[1]


# Changing the pig ids and experimental groups to factors so the ANOVA is done right

# the 'blocking' variable
michalis.constants$Pig                 <- factor(michalis.constants$Pig)

# the independant variable
michalis.constants$Experimental_Groups <- factor(michalis.constants$Experimental_Groups)

# The Vmax values given by the python script are all very small so we need to 
# adjust the units a bit changing Vmax from (mmol/million cells/hr) 
# to (nmol/million cells/hr) this brings the reported values into the same
# units as those commonly reported in other papers
michalis.constants$Vmax <- michalis.constants$Vmax*10^6 

# Using the summarySE found at the top of the script to get the SEs
# and CIs so they can be used in the plotting of errorbars later
michalis.summary.km    <- summarySE(michalis.constants, measurevar="Km", 
                                    groupvars="Experimental_Groups")
michalis.summary.vmax  <- summarySE(michalis.constants, measurevar="Vmax", 
                                    groupvars="Experimental_Groups")
michalis.summary.km2   <- summarySE(michalis.constants, measurevar="Km")
michalis.summary.vmax2 <- summarySE(michalis.constants, measurevar="Vmax")
michalis.summary.r.sq  <- summarySE(michalis.constants, measurevar="R_sq")

# Running the ANOVA to see if there are any signifcant differences
km.anova.values   <- aov(Km ~ Experimental_Groups + Pig, 
                         data = michalis.constants)
vmax.anova.values <- aov(Vmax ~ Experimental_Groups + Pig, 
                         data = michalis.constants)
summary(km.anova.values)
summary(vmax.anova.values)

###############################################################################
# Graphng and displaying the results

#Generating plots for the Km values
km.graph <- ggplot(michalis.summary.km, 
                   aes(Experimental_Groups, 
                       Km, 
                       fill=Experimental_Groups)) + 
  xlab('Oxygen Tension') +
  ylab(bquote(K[m]*' (mM)')) + 
  geom_errorbar(aes(ymin=Km-se, ymax=Km+se),
                width = .2) +
  geom_col(color='black') +
  # ggtitle('Km Vs \nCultured Oxygen Tension') +
  scale_x_discrete(limits = c('2.5%-NP', '5%-NP', '21%-NP'),
                   labels = c('2.5%', '5%', '21%')) +
  scale_fill_manual(values=c('white', 'black', 'grey'))+
  theme_bw() +
  theme(legend.position='none')
plot(km.graph)

# Generating plots for the Vmax values
vmax.graph <- ggplot(michalis.summary.vmax,
                     aes(Experimental_Groups, 
                         Vmax, 
                         fill=Experimental_Groups)) + 
  xlab('Oxygen Tension') + 
  ylab(bquote(V[max]*' (nmol/10'^6*' cells/hr)')) +
  scale_x_discrete(limits = c('2.5%-NP', '5%-NP', '21%-NP'),
                   labels = c('2.5%', '5%', '21%')) +
  geom_errorbar(aes(ymin=Vmax-se, ymax=Vmax+se),
                width = .2) +
  scale_fill_manual(values=c('white', 'black', 'grey'))+
  geom_col(color='black') +
  # ggtitle('Vmax Vs \nCultured Oxygen Tension') +
  theme_bw() +
  theme(legend.position='none')

plot(vmax.graph)
# plot_grid(km.graph, vmax.graph, labels = c('A', 'B'), ncol = 2, nrow = 1)

# Making kinteic plots of Glucose consumption rate vs Glucose concentration.
# The plot will have three lines, one for each experimental condition.
# Mapping the color of each line

# Creating a 'fit' line, the intial concentration needs to be 
# fiddled with to line up with the graph
# This makes sense since it is a constant of integration
sample.curvefit.line <- NULL
for(i in seq(max(sample.curvefit$Glucose.Concentration)*1.2, 0, -0.001)){
  time <- MichalisMenten(i, 
                         sample.curvefit$Glucose.Concentration[1]*1.2, 
                         curve.km, 
                         curve.vmax)
  sample.curvefit.line <- rbind(sample.curvefit.line, c(i,time))
}

#Making the 'fit' datapoints have the same column names
sample.curvefit.line <- data.frame(sample.curvefit.line)
names(sample.curvefit.line) <- c('Glucose.Concentration','Hour')
curve.map <- c('Experimental Data'= 0,
               'Curvefit'=1)

#Plotting the sample values vs the fit line
curvefit.graph <- ggplot(sample.curvefit) +
  geom_point(aes(Hour, 
                 Glucose.Concentration, 
                 shape='Experimental Data'), 
             size=3) +
  expand_limits(y = 5) +
  geom_line(data = sample.curvefit.line, 
            aes(Hour, 
                Glucose.Concentration, 
                linetype='Curvefit')) +
  xlab('Time (Hr)') + 
  ylab('Glucose Concentration (mM)') +
  scale_x_continuous(breaks = seq(0, 12, by = 1)) +
  scale_alpha_manual(values = curve.map) +
  # theme_bw()+
  theme(legend.position='none',
        panel.background = element_rect(fill='grey',
                                        color='black', 
                                        linetype='solid', 
                                        size=1))

# annotate('text', x=0.5, y=0.5, 
#          label='   Km=0.179mM \n
#          Vmax=101(nmol/million cells/hr) 
#          \n R_sq=0.998')
# ggtitle('Curvefit')
plot(curvefit.graph)

# assigning a color to each line
color.map <- c('2.5%' = 'blue', 
               '5%' = 'green', 
               '21%' = 'red',
               'Avg' = 'black') 

# Mapping the linetype of each line
label.map <- c('2.5%' = 'solid',
               '5%' = 'solid', 
               '21%' = 'solid',
               'Avg' = 'solid')

thicc.map <- c('2.5%' = 1,
               '5%' = 1, 
               '21%' = 1,
               'Avg' = 1)
# legend order
label.order <- c('2.5%', '5%', '21%', 'Avg') 

#Generating the plots
kinetic.plots <- ggplot(data.frame(x = c(0, 6), 
                                   y = c(0, max(michalis.summary.vmax$Vmax))), 
                        aes(x, y)) +
  stat_function(fun = function(gluc) kinetics(gluc, 
                                              michalis.summary.vmax, 
                                              michalis.summary.km, 
                                              1),
                aes(color = '2.5%',
                    linetype='2.5%', 
                    size='2.5%')) +
  stat_function(fun = function(gluc) kinetics(gluc, 
                                              michalis.summary.vmax,
                                              michalis.summary.km,
                                              2),
                aes(color = '5%', 
                    linetype='5%', 
                    size='5%')) +
  stat_function(fun = function(gluc) kinetics(gluc,
                                              michalis.summary.vmax,
                                              michalis.summary.km,
                                              3),
                aes(color = '21%',
                    linetype='21%', 
                    size='21%')) +
  stat_function(fun = function (gluc) kinetics(gluc,
                                               michalis.summary.vmax2,
                                               michalis.summary.km2,
                                               1),
                aes(color = 'Avg', 
                    linetype='Avg',
                    size='Avg')) +
  scale_colour_manual('Oxygen\nTension',  
                      values=color.map,
                      breaks=label.order) +
  scale_linetype_manual('Oxygen\nTension', 
                        values=label.map,
                        breaks=label.order) +
  scale_size_manual('Oxygen\nTension', 
                    values=thicc.map, 
                    breaks=label.order) +
  xlab('Glucose Concentration (mM)') + 
  ylab('Glucose Consumption Rate\n(nmol/millioncells/hr)') 
# ggtitle('Michalis-Menten Kinetic Plots')

plot(kinetic.plots)

plot_grid(curvefit.graph, kinetic.plots, 
          labels = c('a', 'b'), 
          ncol = 2, nrow = 1)
# Need to run the lactate_stats.R first
plot_grid(km.graph, vmax.graph, lactate.graph, 
          labels = c('a', 'b', 'c'), 
          ncol = 3, nrow = 1)
