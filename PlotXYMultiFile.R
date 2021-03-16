library("spatstat")
library("matlab")
library("optparse")
library("RJSONIO")
library("dplyr")
library("stringr")
library("ggplot2")
source("/Users/morganoneka/Box/My Stuff/Gcross/functions.R")

# option parsing
option_list = list(
  make_option(c("-r", "--radiusmax"), type="numeric", default=300,
              help="Maximum radius for which to calculate radius [default = %default]"), 
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
  make_option(c("-v", "--verbose"), default=FALSE, help="Verbose mode: phenotype = colname + value"),
  make_option(c("-f", "--filetype"), default="csv", help="File type for input files (default csv)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

radiusmax = opt$radiusmax
xcol = opt$xcol
ycol = opt$ycol
input = opt$input
output = opt$output
jsonphenotype = opt$jsonphenotype
sep = opt$sep
summaryfile = opt$summary
verbose = opt$verbose
filetype = opt$filetype

phenotypes <- fromJSON(jsonphenotype)

for (inx in phenotypes){
  
  for (file in list.files(path=input, pattern=paste("*", filetype, sep="."), full.names=TRUE, recursive=FALSE)){
    print(file)
    XYdata <- read.table(file, header = TRUE, sep = sep, fill = TRUE, row.names=NULL)
    
    ref = inx[[1]][[1]]$Reference
    nonref = inx[[1]][[1]]$`Non-Reference`
    
    ref_phenotypes=getPhenotypeString(ref, verbose)
    nonref_phenotypes=getPhenotypeString(nonref, verbose)
    
    # identify reference and non-reference cells
    XYdata$PastedPhenotype = ""
    try1 <- tryCatch((XYdata[getRowsWithPhenotype(ref,XYdata,verbose), "PastedPhenotype"] = "Reference"), error=function(e) e)
    try2 <- tryCatch((XYdata[getRowsWithPhenotype(nonref,XYdata,verbose), "PastedPhenotype"] = "Non-Reference"), error=function(e) e)
    
    if (inherits(try1, "error") | inherits(try2, "error")){
      next()
    }
    
    
    plot_data = XYdata[which(XYdata$PastedPhenotype %in% c("Reference", "Non-Reference")),]
    plot_data[which(plot_data$PastedPhenotype == "Reference"),"PastedPhenotype"] = gsub("_", " ", ref_phenotypes)
    plot_data[which(plot_data$PastedPhenotype == "Non-Reference"),"PastedPhenotype"] = gsub("_", " ", nonref_phenotypes)
    
    sample <- strsplit(tail(strsplit(file, "/")[[1]], 1), "\\.")[[1]][1]
    
    ggplot(plot_data, aes_string(x=xcol, y=ycol, color="PastedPhenotype")) + geom_point() + xlab("X Position") + ylab("Y Position") + coord_fixed(ratio=1) + ggtitle(paste("Spatial Arrangement for", sample))
    
    filename=paste(paste("sample", sample, paste(ref_phenotypes, nonref_phenotypes, sep="_VS_"), sep="_"), ".png", sep="")
    ggsave(filename, plot = last_plot(), device="png", path=output)
    
  }
  
}

summary_stats = data.frame(matrix(ncol=0, nrow=0))
write.table(summary_stats, summaryfile, sep=",", col.names=TRUE)