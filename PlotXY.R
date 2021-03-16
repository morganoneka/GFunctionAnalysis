library("spatstat")
library("matlab")
library("optparse")
library("RJSONIO")
library("dplyr")
library("stringr")
library("ggplot2")
source("/Users/morganoneka/Box/My Stuff/Gcross/functions.R")

option_list = list(
  make_option(c("-x", "--xcol"), type="character", default="X.Position",
              help="The name of the column giving x position"),
  make_option(c("-y", "--ycol"), type="character", default="Y.Position",
              help="The name of the column giving y position"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input file or directory"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory"),
  make_option(c("-j", "--jsonphenotype"), type="character", default=NULL,
              help="Json file containing cell phenotypes"),
  make_option(c("--sep"), type="character", default=",",
              help="Separator used in input data"),
  make_option(c("--summary"), type="character", default=NULL,
              help="Name for file to save gcross summary info to"),
  make_option(c("-v", "--verbose"), default=FALSE, help="Verbose mode: phenotype = colname + value")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# print(opt)

radiusmax = opt$radiusmax
xcol = opt$xcol
ycol = opt$ycol
input = opt$input
output = opt$output
jsonphenotype = opt$jsonphenotype
sep = opt$sep
verbose = opt$verbose
summaryfile = opt$summary

XYdata <- read.table(input, header = TRUE, sep = sep, fill = TRUE, row.names=NULL)
phenotypes <- fromJSON(jsonphenotype)

if (samplecol != patientcol & ! is.na(samplecol)){
  XYdata[,samplecol] = paste(XYdata[,patientcol], XYdata[,samplecol], sep="_")
}

for (inx in phenotypes){
  XYdata$PastedPhenotype = ""
  
  ref = inx[[1]][[1]]$Reference
  nonref = inx[[1]][[1]]$`Non-Reference`
  
  ref_phenotypes=getPhenotypeString(ref, verbose)
  nonref_phenotypes=getPhenotypeString(nonref, verbose)
  
  XYdata[getRowsWithPhenotype(ref,XYdata,verbose), "PastedPhenotype"] = "Reference"
  XYdata[getRowsWithPhenotype(nonref,XYdata,verbose), "PastedPhenotype"] = "Non-Reference"
  
  
  for (sample in unique(XYdata[,samplecol])){
    sample_data = XYdata[which(XYdata[,samplecol] == sample),]
    plot_data = sample_data[which(sample_data$PastedPhenotype %in% c("Reference", "Non-Reference")),]
    plot_data[which(plot_data$PastedPhenotype == "Reference"),"PastedPhenotype"] = gsub("_", " ", ref_phenotypes)
    plot_data[which(plot_data$PastedPhenotype == "Non-Reference"),"PastedPhenotype"] = gsub("_", " ", nonref_phenotypes)
    ggplot(plot_data, aes_string(x=xcol, y=ycol, color="PastedPhenotype")) + geom_point() + xlab("X Position") + ylab("Y Position") + coord_fixed(ratio=1) + ggtitle(paste("Spatial Arrangement for", sample))
    
    filename=paste(paste("sample", sample, paste(ref_phenotypes, nonref_phenotypes, sep="_VS_"), sep="_"), ".png", sep="")
    ggsave(filename, plot = last_plot(), device="png", path=output)
  }
}

summary_stats = data.frame(matrix(ncol=0, nrow=0))
write.table(summary_stats, summaryfile, sep=",", col.names=TRUE)