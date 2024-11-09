library(IsoformSwitchAnalyzeR)
library(ggplot2)

salmonQuant <- importIsoformExpression(
  parentDir = "path/to/parentDir/where/salmon_outputs/are"  
)

myDesign <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)

#Anno and Fasta have to be gencode v46. 
aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation ="path/to/annotation/gtf_file",   
  isoformNtFasta       = "path/to/transcriptome/fasta_file", 
  fixStringTieAnnotationProblem = F,
  showProgress = T
)

sampleID = colnames(salmonQuant$abundance)[-1]

list2 <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = aSwitchList,
  pathToOutput = paste0('C:/ISOSAR_files/',sampleID,'/'),
  outputSequences      = F, 
  prepareForWebServers = F,
)

#just get a list of top genes

switchList <- extractTopSwitches(
  list2, 
  filterForConsequences = F, #always F
  n = 251,     #Change here to get a longer list for switchNames
  sortByQvals = T
)

print(switchList)

switchNames <- as.data.frame(switchList$gene_name)

print(switchNames) #Top n switches ordered by dIF (diff Isoform Fraction)

#Do this for every gene name in switchNames and compare to literature
switchPlot(
  list2,
  gene='DPP4',   #Pick a gene name from switchNames
  condition1 = 'cnt',
  condition2 = 'IFNa'
)

