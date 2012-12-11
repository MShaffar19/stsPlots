# stsPlots.2.0.R
#
# Author: Meredith Ashby, Lawrence Lee
# Date: 5/12
# Revised 04/12 to accomodate changes in R V2.14.2 and higher, ggplot2_0.9.0, and the replacement of reshape with reshape2
# Revised 12/12 to more robustly handle low quality SMRTcell data
# Description: These are a set of plots generated from the sts.csv file that help in assessing loading of the SMRTcell and pre-mapping quality.

library(ggplot2)
library(methods)
library(reshape2)
library(plyr)


generateStsPlots <- function(stsCsvPath, pdfOutPath) {


  pdf(pdfOutPath, width=11, height=8.5, onefile=T)
  r <- read.table(stsCsvPath, header=TRUE, sep=",", strip.white=T)
  rn <- subset(r, ZmwType=="SEQUENCING")
  rn$HQReadLength <- rn$HQRegionEnd - rn$HQRegionStart

  
  # PAGE 1:  ReadScore heatmap 
  rn$BinnedReadScore <- round_any(rn$ReadScore, 0.75, floor)

  instructions <- "Confirm an even distribution\nof pass / fail Readscores."
  annotation <- paste("ZMWs passing filter:\n", sum(rn$ReadScore >= 0.75), ' ZMWs \n',format(sum(rn$ReadScore >= 0.75) / length(rn$ReadScore) * 100, digits=4), " %\n", sep="")

  
  d <- ggplot(rn, aes(x=X, y=Y)) + geom_tile(aes(fill=factor(BinnedReadScore))) +
      scale_fill_manual(values=c("red", "yellow"), name="Binned Readscore") +
      opts(title="ReadScore Heatmap : Loading / Alignment Metric") +
      coord_cartesian(xlim=c(-175,250)) +
      annotate(geom='text', x=240, y=-170, vjust=0, hjust=1,label=annotation) +  
      annotate(geom='text', x=240, y=240, vjust=1, hjust=1,label=instructions)   

  show(d)


  # PAGE 2:  Productivity heatmap 
  prod0 <- sum(rn$Productivity==0)
  prod1 <- sum(rn$Productivity==1)
  prod2 <- sum(rn$Productivity==2)

  prod0Pct <- paste("Prod=0 ZMWs: ", prod0, " (",format(prod0 / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  prod1Pct <- paste("Prod=1 ZMWs: ", prod1, " (",format(prod1 / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  prod2Pct <- paste("Prod=2 ZMWs: ", prod2, " (",format(prod2 / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  activePct <- paste("Total Active ZMWs: ",(prod1+prod2), " (", format((prod1 + prod2) / length(rn$Productivity) * 100, digits=4), "%)", sep="")
  instructions <- "Poisson distribution indicates that\nProd=1 ZMWs are maximized at 36.7%\nif 63.2% of total ZMWs are active.\n"
  annotation <- paste(prod0Pct, '\n', prod1Pct, '\n', prod2Pct, '\n', activePct, sep="")
  
  d <- ggplot(rn, aes(x=X, y=Y)) + geom_tile(aes(fill=factor(Productivity))) +
      scale_fill_manual(values=c("red","yellow","blue"),name="Productivity") +
      opts(title="Productivity Heatmap : Loading Metric") +
      coord_cartesian(xlim=c(-175,250)) +
      annotate(geom='text', x=240, y=min(rn$Y), vjust=0, hjust=1,label=annotation) +
      annotate(geom='text', x=240, y=max(rn$Y), vjust=1, hjust=1,label=instructions)
  
  show(d)

  
  # PAGE 3: ReadScore Distribution vs. Productivity
  
  d <- ggplot(subset(rn, rn$ReadScore > 0.1), aes(ReadScore, fill=factor(Productivity))) +
         geom_bar() +
         scale_fill_hue(name="Productivity") +
         opts(title="ReadScore Distribution by Productivity : Sequencing Quality Metric") +
	 geom_vline(xintercept = 0.75, color="red") +
         facet_grid(Productivity~.)  

  show(d)


  # PAGE 4: HQReadLength boxplot 
  p <- ggplot(rn[rn$ReadScore > 0.1,], aes(x=factor(round_any(HQReadLength, 500)), y=ReadScore, color=as.factor(Productivity))) +
       opts(axis.text.x=theme_text(angle=-90, hjust=0), title="ReadScore Distribution by HQ Readlength : Sequencing Quality Metric") +
       geom_boxplot(outlier.size=0) + xlab("Binned HQ ReadLength") +
       scale_colour_hue(name="Productivity") 
  show(p)

  
  # PAGE 5:Comparison of raw vs HQ region only readlength distributions - How much of the sequencing-zmw data is crap?
  hq <- melt(rn, measure=c("NumBases","HQReadLength"), id=c("X","Y","ZmwType","ReadScore"))
  colnames(hq) <- c("X","Y","ZmwType","ReadScore","Raw_vs_HQ_RL","value")
  instructions <- "A large discrepancy between raw and\nHQ read lengths indicates noisy ZMWs.\n\nCompare HQ read length plot to expected\n readlength values for your movie length."

  p <- qplot(data=subset(hq, hq$ReadScore > 0.1), x=value, geom='freqpoly', color=Raw_vs_HQ_RL ) +
    xlab("Length in Bases") +
    ylab("Counts") +
    opts(title="Raw vs HQ Region Readlength Distributions : Sequencing Quality Metric") +
    scale_colour_hue(name="Raw vs. HQ Readlength") +
    annotate(geom='text', x=max(hq$value), y=0, vjust=0, hjust=1,label=instructions)

  show(p)

  
  # PAGE 6: What fraction of called bases are in an HQ region?
  rn$HQFraction <- rn$HQReadLength / rn$NumBases
  instructions <- "A large discrepancy between raw and\nHQ readlengths indicates noisy ZMWs."
  p <- qplot(data=subset(rn, ReadScore > 0.1), x=HQFraction, geom='histogram') +
    xlab("Fraction of Called Bases that are in an HQRegion") +
    ylab("Counts") +
    opts(title="Per ZMW Sequencing / Noise Ratio") +
    annotate(geom='text', x=0, y=500, vjust=0, hjust=0, label=instructions)

  show(p)

  
  # PAGE 7: SNR heatmap
  mb <- melt(rn, measure.var=c("SnrMean_T", "SnrMean_G", "SnrMean_A", "SnrMean_C"), id.var=c("X", "Y", "ZmwType", "ReadScore"))
  q1 <- quantile(subset(mb, ZmwType=="SEQUENCING")$value, 0.01, na.rm=TRUE)
  q99 <- quantile(mb$value, 0.99, na.rm=TRUE)
  mb$Baseline <- pmax(0.9*q1, pmin(mb$value, q99))

  p <- qplot(X,Y, fill=Baseline, geom="tile", data=mb, main = "Mean SNR Heatmap : Alignment Metric") +
    facet_wrap(~ variable) +
    scale_fill_gradientn(colours = rainbow(7))
  
  show(p)
  

  # PAGE 8: SNR distributions
  #instructions <- "MINIMUM: Channel T > 4, Channel A > 5.5\nOPTIMUM: Channel T 4.5-7, Channel A 7-10\n(FCR, IMT chips)"
  instructions <- data.frame(variable=c('SnrMean_T','SnrMean_G','SnrMean_A','SnrMean_C'), Annotation=c('Minimum: 4.0\nOptimum: 4.5-7', '','Minimum: 5.5\nOptimum: 7-10',''))

  p <- qplot(data=subset(mb, ReadScore > 0.1), x=value, geom='freqpoly', color=variable) +
    xlab("ZMW Mean SNR") + ylab("Counts") +
    opts(title="Per Channel Mean SNR Distribution")  +
    facet_wrap(~variable) + coord_cartesian(xlim=c(0,16)) +
    scale_x_continuous(breaks=seq(2,14,2)) +
    scale_colour_hue(name="Channel SNR") +
    geom_text(aes(label=Annotation), data=instructions, x=1, y=0, vjust=0, hjust=0, show_guide=FALSE)
  
  show(p)
  
  
  # PAGE 9: IPD histogram : oxygen exclusion test

  prod <- subset(rn,ReadScore >= 0.75 & Productivity == 1 & BaseIpd < 2)
  ipdMean <- mean(prod$BaseIpd)
  interpretation <- ifelse(ipdMean<0.25,'Normal','Slow')
  annotation <- paste('Mean IPD(<0.25 is normal): ', round_any(ipdMean,0.01), 's\nnReads(>22k is normal): ',dim(prod)[1], '\nMode of distribution should\n be to the left of vertical line\n\n',interpretation, '\n',sep="")
  
  p <- qplot(data=prod,x=BaseIpd,y=..count..,binwidth=0.02,geom='freqpoly',main='IPD Histogram : Oxygen Exclusion') +
    theme_bw() +
    geom_vline(xintercept=0.25,lwd=1,color=ifelse(ipdMean<0.25,'blue','red')) +
    scale_x_continuous("Base IPD (s)") +
    annotate(geom='text', x=max(prod$BaseIpd), y=0, vjust=0, hjust=1,label=annotation)
  
  show(p)

  dev.off()
}

plotStsCsvFolder <- function(stsCsvFolder) {
  outputFolder <- file.path(stsCsvFolder, "Analysis_Results")
  tryCatch(dir.create(outputFolder), simpleWarning = function(e) stop(paste("Unable to create directory: ", outputFolder)))

  csvFiles <- list.files(path = stsCsvFolder, pattern = "*.sts.csv")
  print(csvFiles)

  for (stsCsvFile in csvFiles) {
    print(stsCsvFile)
    pdfOutput <- file.path(outputFolder, paste(strsplit(stsCsvFile, split='\\.csv')[[1]][1], '.pdf', sep=""))
    generateStsPlots(file.path(stsCsvFolder,stsCsvFile), pdfOutput)
  }
}


runStsPlots <- function(args) {

  if (args[1] !="-file" & args[1]!="-folder")
    print("stsPlots requires either the -file or -folder flag  Ex: Rscript --vanilla stsPlots.2.0.R -file yourfile.sts.csv")

  inputType <- args[1]

  if (inputType=="-file") {
    if (length(args)!=2) print("You must supply a file name as a second argument.")
    stsCsvFile <- args[2]
    pdfOutput <- paste(strsplit(stsCsvFile, split='\\.csv')[[1]][1],'.pdf',sep="")
    generateStsPlots(stsCsvFile,pdfOutput)
  }

  if (inputType == "-folder") {
    if (length(args)!=2) print("You must supply a folder path as a second argument")
    stsCsvFolder <- args[2]
    plotStsCsvFolder(stsCsvFolder)
  }
}

args <- commandArgs(trailingOnly = T)
runStsPlots(args)

