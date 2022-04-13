# This script is for finding the differences between oxygen tensions for NP
# cell gene expression after a 24 hour adaptation period. A correlational 
# analysis between the transcription factor HIF-1a and all other genes of 
# interest was also carried out.

# Author <- Lukas Jaworski, l.jaworski@umiami.edu
# Institution <- University of Miami
# PI <- Alicia Jackson, a.jackson2@miami.edu

# Copyright University of Miami 2017
###############################################################################
## Library imports
###############################################################################

library(ReadqPCR)
library(NormqPCR)
library(ggplot2)
library(gridExtra)
library(cowplot)

###############################################################################
## Functions
###############################################################################

aov.results <- function(genes, data, formula){
  # A function to get the statistically significant results
  anova.results <- list()
  for (x in genes){
    anova.results[[x]] <- aov(as.formula(gsub('{}', x, formula, fixed=T)),
                              data = data)
  }
  
  # ANOVA results and perform post hoc tests
  aov.summary<-list()
  for (x in genes){
    temp.summary <- summary.aov(anova.results[[x]])[[1]]
    aov.summary[[x]][[1]] <- temp.summary[temp.summary['Pr(>F)']<=0.05 &
                                          !is.na(temp.summary['Pr(>F)']),
                                          'Pr(>F)',
                                          drop=F]
    if (nrow(aov.summary[[x]][[1]])>0){
      temp.tukey <- TukeyHSD(anova.results[[x]])
      aov.summary[[x]][[2]] <- lapply(temp.tukey,
                                      function(x) x[x[,'p adj']<=0.05,'p adj'])
    }
  }
  return(aov.summary)
}

###############################################################################
## Script
###############################################################################
# Retrieving and arranging the data

# Set the working directory to where the data is kept, make it workstation 
# agnostic
if (Sys.info()[['sysname']]=='Linux'){
  working.directory <- paste0('/home/',Sys.info()[['user']],
                             '/Dropbox/Lab/Experiments/PCR/',
                             'a3-gcr(ox)-NP.AF/Script')
} else {
working.directory <- paste0('C:/Users/',Sys.info()[['user']],
                            '/Dropbox/Lab/Experiments/PCR/',
                            'a3-gcr(ox)-NP.AF/Script')
}
working.directory <- paste0('P:/Dropbox/Lab/Experiments/PCR/',
                            'a3-gcr(ox)-NP.AF/Script')
setwd(working.directory)

# Read in the Cq values and run the deltaCq method on them
qPCRBatch.Cq <- read.qPCR('pcrfile.txt')
qPCRBatch.dCq <- deltaCq(qPCRBatch = qPCRBatch.Cq, hkgs = "GAPDH")

# Getting the indices of Oxygen levels, set a column with O2 levels as factors
# This is our independant variable
O2groups <- list()
O2groups[[1]] <- grep('_2.5_', sampleNames(qPCRBatch.dCq))
O2groups[[2]] <- grep('_21_', sampleNames(qPCRBatch.dCq))
O2groups[[3]] <- grep('_5_', sampleNames(qPCRBatch.dCq))
names(O2groups) <- c('low', 'high', 'phys')

# Getting the indices of Pig source, set a column with pigs as factors
# This is our blocking variable
pigs <- list()
pigs[[1]]<-grep('P19',sampleNames(qPCRBatch.dCq))
pigs[[2]]<-grep('P20',sampleNames(qPCRBatch.dCq))
pigs[[3]]<-grep('P21',sampleNames(qPCRBatch.dCq))
pigs[[4]]<-grep('P22',sampleNames(qPCRBatch.dCq))
pigs[[5]]<-grep('P23',sampleNames(qPCRBatch.dCq))
names(pigs) <- c('p19','p20','p21','p22','p23')

gene.expression<-data.frame(t(exprs(qPCRBatch.dCq)))
gene.expression$GAPDH <- NULL
#This will cycle through all comparisons low-high, low-phys, and high-low
gene.expression$level <- NA
gene.expression$level[O2groups$low]<-'low'
gene.expression$level[O2groups$phys]<-'phys'
gene.expression$level[O2groups$high]<-'high'
gene.expression$level<-as.factor(gene.expression$level)

gene.expression$pig <- NA
gene.expression$pig[pigs$p19]<-'p19'
gene.expression$pig[pigs$p20]<-'p20'
gene.expression$pig[pigs$p21]<-'p21'
gene.expression$pig[pigs$p22]<-'p22'
gene.expression$pig[pigs$p23]<-'p23'
gene.expression$pig <- as.factor(gene.expression$pig)
# kruskal.test(GLUT.1 ~ level, gene.expression)

name.table <- list(Agg='Aggrecan', 
                   Col_IA1='Collagen I', 
                   Col_IIA1='Collagen II',
                   GAPDH='GAPDH',
                   MMP_1='MMP 1',
                   MMP_13='MMP 13',
                   MMP_3='MMP 3',
                   MT_1='MT 1',
                   'T'='Brachury',
                   TIMP_1='TIMP 1',
                   TIMP_2='TIMP 2',
                   TIMP_3='TIMP 3',
                   GLUT.1='GLUT-1',
                   GLUT.3='GLUT-3',
                   GLUT.9='GLUT-9',
                   HIF.1=expression(paste('HIF 1',alpha)))

gene.names <- setdiff(names(name.table),'GAPDH')

# run the anovas
anova.summary.o2 <- aov.results(gene.names,
                                gene.expression,
                                '{} ~ level + pig')
# pvals <- data.frame(row.names = featureNames(qPCRBatch.dCq))
# for (x in 1:(length(names(O2groups))-1)){
#   for (y in (x+1):length(names(O2groups))){
#     for(i in row.names(pvals)) {
#       case <- exprs(qPCRBatch.dCq)[i,O2groups[[x]]]
#       control <- exprs(qPCRBatch.dCq)[i,O2groups[[y]]]
#       if(sum(is.na(case)) != 0 | sum(is.na(control)) != 0) { 
#         # gets rid of NaN values (mostly the housekeeping gene)
#         pvals[i, paste(names(O2groups[x]),
#                        names(O2groups[y]),
#                        sep = '_')] <- NA
#       }else{
#         pvals[i, paste(names(O2groups[x]),
#                        names(O2groups[y]),
#                        sep = '_')] <- wilcox.test(case,control)$p.value
#       }
#     }
#   }
# }
# # nested for loops, does all the comparisons I want for the exploratory 
# # data analysis.
# # (Looking to see if there is consistent trend with each pig or not)
# 
# for (z in 1:length(names(pigs))){
#   for (x in 1:(length(names(O2groups))-1)){
#     for (y in (x+1):length(names(O2groups))){
#       for(i in row.names(pvals)) {
#         case <- exprs(qPCRBatch.dCq)[i,intersect(pigs[[z]],O2groups[[x]])]
#         control <- exprs(qPCRBatch.dCq)[i,intersect(pigs[[z]],O2groups[[y]])]
# 
#         # gets rid of NaN values (mostly the housekeeping gene)
#         if(sum(is.na(case)) != 0 | sum(is.na(control)) != 0) {
#           pvals[i, paste(names(pigs[z]),
#                          names(O2groups[x]),
#                          names(O2groups[y]),
#                          sep = '_')] <- NA
#         }else{
#           pvals[i, paste(names(pigs[z]),
#                          names(O2groups[x]),
#                          names(O2groups[y]),
#                          sep = '_')] <- wilcox.test(case,control)$p.value
#         }
#       }
#     }
#   }
# }

#Adjusting p-values for multiple comparisons
pvals_to_adjust<-rbind(t(t(pvals[[1]])), 
                       t(t(pvals[[2]])), 
                       t(t(pvals[[3]])))
adjPVals <- p.adjust(as.numeric(pvals_to_adjust[,1]), method = "BY")
# show old p-vals next to adjust values
pvals_to_adjust <- cbind(pvals_to_adjust,adjPVals) 
# It adjusted them all to 1 so...
colSums(pvals < 0.05, na.rm = TRUE)
# Running the ANOVA for sensitivity
anova.summary <- aov.results(gene.names,gene.expression,'{}~level')

#finding correlations between HIF transcription factor and other genes studied
dCq <- exprs(qPCRBatch.dCq)

dCq.t <- t(dCq)
dCq.t <- dCq.t[, -4]
p.correlations<-NULL
for (x in 1:15){
  p.correlations <- rbind(p.correlations, 
                          c(colnames(dCq.t)[x], 
                            cor.test(dCq.t[,7], 
                                     dCq.t[,x], 
                                     method = 'pearson')$estimate, 
                            cor.test(dCq.t[,7],
                                     dCq.t[,x], 
                                     method = 'pearson')$p.value))
}
p.correlations <- p.correlations[-7, ]
adjPVals.correl <- p.adjust(as.numeric(p.correlations[,3]),
                            method = "bonferroni")
p.correlations <- cbind(p.correlations, adjPVals.correl)

#graphing dCqs and correlation plots
#plotting corraltional values

t.hif.graph <- ggplot(gene.expression, aes(HIF.1,T)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('Brachyury')

col1.hif.graph <- ggplot(gene.expression, aes(HIF.1,Col_IA1)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('Collagen I')

col2.hif.graph <- ggplot(gene.expression, aes(HIF.1,Col_IIA1)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('Collagen II')

agg.hif.graph <- ggplot(gene.expression, aes(HIF.1,Agg)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('Aggrecan')

anabolic.hif.graphs <- plot_grid(col1.hif.graph, 
                                 col2.hif.graph, 
                                 agg.hif.graph, 
                                 ncol = 3,
                                 labels = c('i','ii','iii'),
                                 hjust = -0.6)

glut1.hif.graph <- ggplot(gene.expression, aes(HIF.1,GLUT.1)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('GLUT 1')

glut3.hif.graph <- ggplot(gene.expression, aes(HIF.1,GLUT.3)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('GLUT 3')

glut9.hif.graph <- ggplot(gene.expression, aes(HIF.1,GLUT.9)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('GLUT 9')

gluts.hif.graphs <- plot_grid(glut1.hif.graph,
                              glut3.hif.graph, 
                              glut9.hif.graph,
                              ncol = 3,
                              labels = c('i','ii','iii'),
                              hjust = -0.6)

mt1.hif.graph <- ggplot(gene.expression, aes(HIF.1,MT_1)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('MT 1')

mmp1.hif.graph <- ggplot(gene.expression, aes(HIF.1,MMP_1)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('MMP 1')

mmp3.hif.graph <- ggplot(gene.expression, aes(HIF.1,MMP_3)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('MMP 3')

mmp13.hif.graph <- ggplot(gene.expression, aes(HIF.1,MMP_13)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('MMP 13')

catabolic.hif.graphs <- plot_grid(mt1.hif.graph,
                                  mmp1.hif.graph,
                                  mmp3.hif.graph, 
                                  mmp13.hif.graph,
                                  ncol = 4,
                                  labels = c('i','ii','iii', 'iv'),
                                  hjust = -0.6)

timp1.hif.graph <- ggplot(gene.expression, aes(HIF.1,TIMP_1)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('TIMP 1')
timp2.hif.graph <- ggplot(gene.expression, aes(HIF.1,TIMP_2)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('TIMP 2')
timp3.hif.graph <- ggplot(gene.expression, aes(HIF.1,TIMP_3)) + geom_point() + 
  geom_smooth(method = 'lm',se = FALSE) +
  xlab(expression(paste('HIF 1', alpha))) +
  ylab('TIMP 3')

timps.hif.graphs <- plot_grid(timp1.hif.graph,
                              timp2.hif.graph,
                              timp3.hif.graph,
                              ncol = 3,
                              labels = c('i','ii','iii'),
                              hjust = -0.6)

all.graphs <- plot_grid(t.hif.graph,
                        anabolic.hif.graphs,
                        catabolic.hif.graphs, 
                        timps.hif.graphs, 
                        gluts.hif.graphs,
                        nrow = 5,
                        labels = c('A','B','C','D','E'), 
                        hjust = 0.15, 
                        vjust = 0.1)

#Used to label the different groups
exp_grps <- c(rep('Low', length(O2groups$low)), 
              rep('Phys',length(O2groups$phys)),
              rep('High',length(O2groups$high)))



#The next three functions just make plotting easier

theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.line = element_line(colour = "black"),
             text = element_text(size = 22))

arrangeData <- function(gene){ 
  # isolate the data for a single gene and label its group
  vals<-c(dCq[gene, O2groups$low], 
          dCq[gene, O2groups$phys], 
          dCq[gene, O2groups$high])
  return(vals)
}

generatePlot <- function(df, gene){
  graph <- ggplot(df, aes(Level, -dCqs)) +
    geom_boxplot() +
    scale_x_discrete(limits = c('High','Phys','Low'),
                     labels = c('21%', '5%', '2.5%')) +
    xlab('Oxygen Tension') +
    ylab(expression(paste('-', Delta, 'Cq'))) +
    labs(title = name.table[[gene]])
  return(graph)
}

plotter <- function(gene){
  vals <- arrangeData(gene)
  df <- data.frame(dCqs=vals,Level=exp_grps)
  return(generatePlot(df, gene))
}

#I need to store all the plots so they can be nicely arranged
for (x in gene.names){
  assign(paste0(tolower(x),'.graph'), plotter(x))
}

anabolic.dCq.graphs <- plot_grid(agg.graph,
                                 col_ia1.graph, 
                                 col_iia1.graph,
                                 ncol = 3, 
                                 labels = 'AUTO',
                                 hjust = 0.5)
plot(anabolic.dCq.graphs)

catabolic.dCq.graphs <- plot_grid(mmp_1.graph, 
                                  mmp_3.graph,
                                  mmp_13.graph, 
                                  mt_1.graph,
                                  ncol = 4, 
                                  labels = 'AUTO')
plot(catabolic.dCq.graphs)

timps.dCq.graphs <- plot_grid(timp_1.graph, 
                              timp_2.graph, 
                              timp_3.graph,
                              ncol = 3,
                              labels = 'AUTO')
plot(timps.dCq.graphs)

gluts.dCq.graphs <- plot_grid(glut.1.graph, 
                              glut.3.graph, 
                              glut.9.graph,
                              ncol = 3, 
                              labels = 'AUTO')
plot(gluts.dCq.graphs)

transcript.dCq.graphs <- plot_grid(t.graph,
                                   hif.1.graph,
                                   ncol = 2, 
                                   labels = 'AUTO')
plot(transcript.dCq.graphs)


## To make and save dCq plots for a quick presentation
# for(i in rownames(dCq)){
#   vals<-c(dCq[i, O2groups$low], 
#           dCq[i, O2groups$phys], 
#           dCq[i, O2groups$high])
#   df <- data.frame(dCqs=vals,Level=exp_grps)
#   graph <- ggplot(df, aes(Level, dCqs)) +
#     geom_boxplot() +
#     geom_point() +
#     labs(title = i)
#   png(paste0(i,'.png'))
#   plot(graph)
#   dev.off()
# }
# for(x in 1:length(names(pigs))){
#   exp_grps<-c(rep('Low',length(intersect(O2groups$low,pigs[[x]]))),
#               rep('Phys',length(intersect(O2groups$phys,pigs[[x]]))),
#               rep('High',length(intersect(O2groups$high,pigs[[x]]))))
#   pig<-names(pigs)[x]
#   for(i in rownames(dCq)){
#     vals<-c(dCq[i,intersect(O2groups$low,pigs[[x]])],
#             dCq[i, intersect(O2groups$phys,pigs[[x]])],
#             dCq[i, intersect(O2groups$high,pigs[[x]])])
#     df <- data.frame(dCqs=vals,Level=exp_grps)
#     graph<-ggplot(df, aes(Level,dCqs)) +
#       geom_boxplot() +
#       geom_point() +
#       labs(title = paste(pig, i))
#     png(paste0(pig,i,'.png'))
#     plot(graph)
#     dev.off()
#   }
# }

##Starting the ddCq method to do delta delta Cq

##Setting up contrast matrices(CM)
# ones <- rep(1,length(O2groups$phys))
# physO2.CM <- rep(0,61)
# physO2.CM[O2groups$phys] <- ones
# ones<- rep(1,length(O2groups$high))
# highO2.CM <- rep(0,61)
# highO2.CM[O2groups$high] <- ones
# ones<- rep(1,length(O2groups$low))
# lowO2.CM <- rep(0,61)
# lowO2.CM[O2groups$low] <- ones
# 
# contM <- cbind(physO2.CM,highO2.CM) # make a contrast matrix
# colnames(contM) <- c("control","case")
# rownames(contM) <- sampleNames(qPCRBatch.Cq)
# ddCq.out <- deltaDeltaCq(qPCRBatch = qPCRBatch.Cq, maxNACase = 0, 
#                          maxNAControl = 0, hkgs="GAPDH", contrastM=contM, 
#                          case="case", control = "control")

