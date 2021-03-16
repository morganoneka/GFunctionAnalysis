library("spatstat")
library("matlab")
library("optparse")
library("RJSONIO")
library("dplyr")
library("stringr")
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

# if (sep == "tab"){
#   sep = "\t"
# }


# read in phenotype json file
phenotypes <- fromJSON(jsonphenotype)

summary_stats = data.frame(matrix(ncol=2, nrow=0), stringsAsFactors = FALSE)

# iterate through phenotypes
for (inx in phenotypes){
  
  # get reference and non-reference labels from each entry of phenotype json
  ref = inx[[1]][[1]]$Reference
  nonref = inx[[1]][[1]]$`Non-Reference`
  
  # create identifying strings for each phenotype
  ref_phenotypes=getPhenotypeString(ref, verbose)
  nonref_phenotypes=getPhenotypeString(nonref, verbose)
  
  # create the name of the file for gcross output
  filename = paste(paste(ref_phenotypes, nonref_phenotypes, sep="_VS_"), "_THEORETICAL", ".csv", sep="")
  
  # create data frame to save results
  gcross_results = data.frame(matrix(ncol=radiusmax+2, nrow=0))
  
  # iterate thru each file in the input directory
  for (file in list.files(path=input, pattern=paste("*", filetype, sep="."), full.names=TRUE, recursive=FALSE)){
    print(file)
    
    # read in file
    XYdata <- read.table(file, header = TRUE, sep = sep, fill = TRUE, row.names=NULL)
    XYdata <- XYdata[complete.cases(XYdata),]
    
    # identify reference and non-reference cells
    XYdata$PastedPhenotype = ""
    try1 <- tryCatch((XYdata[getRowsWithPhenotype(ref,XYdata,verbose), "PastedPhenotype"] = "Reference"), error=function(e) e)
    try2 <- tryCatch((XYdata[getRowsWithPhenotype(nonref,XYdata,verbose), "PastedPhenotype"] = "Non-Reference"), error=function(e) e)
    
    if (inherits(try1, "error") | inherits(try2, "error")){
      next()
    }
    
    # save number of reference and non-reference cells to summary stats
    summary_stats <- rbind(summary_stats, c(ref_phenotypes, as.numeric(nrow(XYdata[which(XYdata$PastedPhenotype == "Reference"),]))), stringsAsFactors=FALSE)
    summary_stats <- rbind(summary_stats, c(nonref_phenotypes, as.numeric(nrow(XYdata[which(XYdata$PastedPhenotype == "Non-Reference"),]))), stringsAsFactors=FALSE)
    
    # point process stuff
    cw <- convexhull.xy(XYdata[,xcol], XYdata[,ycol])
    ww <- owin(poly = cw$bdry)
    pp <- as.ppp(cbind(unlist(XYdata[,xcol]),unlist(XYdata[,ycol])), W = c(-1,2*max(XYdata[,xcol]),-1,2*max(XYdata[,ycol])))
    pp <- pp %mark% factor(unlist(XYdata$PastedPhenotype))
    pp$window <- ww
    
    if (sum(unlist(XYdata$PastedPhenotype) == "Reference") < 1 | sum(unlist(XYdata$PastedPhenotype) == "Non-Reference") < 1) {
      gcross_results = rbind.data.frame(gcross_results, c(file, rep(0, radiusmax+1)), stringsAsFactors = FALSE)
    } else{
      Gcross1 = Gcross(pp,"Reference","Non-Reference",r = 0:radiusmax)
      gcross_results = rbind.data.frame(gcross_results, c(file, t(Gcross1$theo)), stringsAsFactors = FALSE)
    }
  }
  
  # write out gcross results for this interaction
  write.table(gcross_results, gsub("//","/",fullfile(output, filename)), sep=",", col.names=FALSE, quote=FALSE)
}

colnames(summary_stats) <- c("CellType", "Count")
summary_stats <- distinct(summary_stats)
write.table(summary_stats, summaryfile, sep=",", col.names=TRUE, quote=FALSE, row.names=FALSE)


