#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#User need to supply the following: the path to the swept files, args[1], the number of cores for parallelization, args[2], and the name of the chromosome (or contig) for which the plot is generated, args[3]
if (!require(gsalib)){
 install.packages('gsalib')
 library(gsalib)
}
if (!require(dplyr)) {
 install.packages('dplyr') 
 library(dplyr)
}
if (!require(doParallel)){
 install.packages('doParallel')
 library('doParallel')
}
if (!require(foreach)){
 install.packages('foreach')
 library(foreach)
}
#if (!require(doSNOW)){
# install.packages('doSNOW')
# library(doSNOW)
#}
# for submitting jobs to a cluster... Seems too complicated for me, so I'm opting for something else!
#setwd("/home/groups/hpcbio_shared/azza/GIAB/src")
#library(BatchJobs)
#loadConfig()

parameters = data.frame(
  symbol=c('bqsrBAQGOP', 'ddq', 'ics', 'idq', 'lqt', 'maxCycle', 'mcs', 'mdq'), 
  meaning=c('BAQ Gap Open Penalty', 'Deletions default quality', 
            'Indels Context Size', 'Insertions default quality', 'Low quality tail',
            'Max cycle value','Mismatches Context Size',
            'Mismatches default quality')
)
defaults = c(40, 45, 3, 45, 2, 500, 2, -1)

setwd(paste0(args[1],"/default")) #go to the default folder
before = Sys.glob(paste0('*.',args[3],'.recal_table.*'))
after = Sys.glob(paste0('*.', args[3], '.after_recal_table.*')) 
db = gsa.read.gatkreport(before)
qb = db$RecalTable1
da = gsa.read.gatkreport(after)
qa = da$RecalTable1
fwse_before.default = sum((qb$EmpiricalQuality - qb$QualityScore)^2*qb$Observations/sum(as.numeric(qb$Observations)),na.rm = T)
fwse_after.default = sum((qa$EmpiricalQuality - qa$QualityScore)^2*qa$Observations/sum(as.numeric(qa$Observations)), na.rm = T)

dirs = list.dirs('../')  
dirs = dirs[-which(dirs=='..//default')]
param = gsub('..//','',dirs)[2:length(dirs)]
param_name = as.character(parameters[parameters[[1]]==param,2])
####################################################################

   #registerDoMC()
cluster = makeCluster(args[2]) #Second argument is number of core to use
   #registerDoSNOW(cluster)
registerDoParallel(cluster)

foreach (i = 1:length(dirs)+1 ) %dopar% {
###################################################################
#for (i in seq_along(dirs)+1 ) {
  library(gsalib)
  library(dplyr)
  pdffile=paste0(args[1],'BQSR_sweep_plot_on_',args[3],'_',param_name[i-1],'.pdf') #args[3] is the chrmosome name
  pdf(pdffile,title = 'Effect of parameter chages in the bqsr result')
  setwd(dirs[i])
  before =Sys.glob(paste0('*.',args[3],'.recal_table.*')) 
  after =Sys.glob(paste0('*.', args[3], '.after_recal_table.*')) 

  values=c()
  fwsebefore=c()
  fwseafter=c()
  for (file in seq_along(before)){
    valueb = gsub('[^0-9]', '', before[file])
    t = paste0('\\b',valueb,'\\b')
    valuea = grep(t,after)
    rm(t)
    
    if (length(valuea) ){
      db = gsa.read.gatkreport(before[file])
      qb = db$RecalTable1  #The bqsr table of the reported quality as covariate
      da = gsa.read.gatkreport(after[valuea])
      qa = da$RecalTable1
      fwse_before = sum((qb$EmpiricalQuality - qb$QualityScore)^2*qb$Observations/sum(as.numeric(qb$Observations)),na.rm = T)
      fwse_after = sum((qa$EmpiricalQuality - qa$QualityScore)^2*qa$Observations/sum(as.numeric(qa$Observations)), na.rm = T)
      values = c(values, valueb)
      fwsebefore = c(fwsebefore, fwse_before)
      fwseafter = c(fwseafter, fwse_after)
    } 
  }
  FWSE = data.frame(values,fwsebefore,fwseafter)
  colnames(FWSE) = c('value','before','after')
  FWSE$value = as.double(as.character(FWSE$value))
  FWSE = arrange(FWSE,value)
  filename=paste0(args[1],"/fwsevalues",param_name[i-1],".txt")
  write.table(FWSE,filename,sep='\t',append=TRUE, col.names=TRUE)
  ## plot fwse changes Vs parameter changes; i.e as you vary each parameter, 
  #show what happens to fwse
  # actually, you need to know the relative improvement in quality scores after bqsr, 
  # which is |fwse_before -  fwse_after|
  par(xpd=T, mar=c(5,4,4,6))
  plot(FWSE$value, abs(FWSE$before-FWSE$after), 
       main = paste('FWSE after in response to \n "',param_name[i-1],'" change'),
       xlab = param[i-1], ylab = '|change in QUAL|')
  lines(FWSE$value, abs(FWSE$before-FWSE$after))
  points(defaults[i-1],abs(fwse_before.default-fwse_after.default), pch = 19, col ='red')
  legend('topleft',c('default'), pch = 19,
         col=c('red'),inset=c(1.03,0), cex = 0.6 )
  dev.off()
}


stopCluster(cluster)
