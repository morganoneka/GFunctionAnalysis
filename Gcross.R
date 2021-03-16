library("spatstat")
library("matlab")
library("optparse")
library("RJSONIO")
library("dplyr")
library("stringr")
# source("/Users/morganoneka/Box/My Stuff/Gcross/functions.R")
source("./functions.R")

# option parsing
option_list = list(
  make_option(c("-r", "--radiusmax"), type="numeric", default=300,
              help="Maximum radius for which to calculate radius [default = %default]"),
  make_option(c("-x", "--xcol"), type="character", default="X.Position",
              help="The name of the column giving x position"),
  make_option(c("-y", "--ycol"), type="character", default="Y.Position",
              help="The name of the column giving y position"),
  make_option(c("--patientcol"), type="character", default="Patient",
              help="The name of the column giving patient number/ID"),
  make_option(c("--samplecol"), type="character", default="Sample",
              help="The name of the column giving sample number/ID"),
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
  make_option(c("-w", "--imagewise"), default=FALSE, help="Use this flag for image-wise gcross")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# print(opt)

radiusmax = opt$radiusmax
xcol = opt$xcol
ycol = opt$ycol
patientcol = opt$patientcol
samplecol = opt$samplecol
input = opt$input
output = opt$output
jsonphenotype = opt$jsonphenotype
sep = opt$sep
summaryfile = opt$summary
verbose = opt$verbose
imagewise = opt$imagewise

XYdata <- read.table(input, header = TRUE, sep = sep, fill = TRUE, row.names=NULL)
phenotypes <- fromJSON(jsonphenotype)

summary_stats = data.frame(matrix(ncol=2, nrow=0), stringsAsFactors = FALSE)

for (inx in phenotypes){
  XYdata$PastedPhenotype = ""

  ref = inx[[1]][[1]]$Reference
  nonref = inx[[1]][[1]]$`Non-Reference`

  ref_phenotypes=getPhenotypeString(ref, verbose)
  nonref_phenotypes=getPhenotypeString(nonref, verbose)

  XYdata[getRowsWithPhenotype(ref,XYdata,verbose), "PastedPhenotype"] = "Reference"
  XYdata[getRowsWithPhenotype(nonref,XYdata,verbose), "PastedPhenotype"] = "Non-Reference"

  summary_stats <- rbind(summary_stats, c(ref_phenotypes, as.numeric(nrow(XYdata[which(XYdata$PastedPhenotype == "Reference"),]))), stringsAsFactors=FALSE)
  summary_stats <- rbind(summary_stats, c(nonref_phenotypes, as.numeric(nrow(XYdata[which(XYdata$PastedPhenotype == "Non-Reference"),]))), stringsAsFactors=FALSE)

  filename=""
  if (imagewise){
    filename=paste(paste(paste(ref_phenotypes, nonref_phenotypes, sep="_VS_"), "imagewise", sep="_"), ".csv", sep="")
  } else{
    filename=paste(paste(paste(ref_phenotypes, nonref_phenotypes, sep="_VS_"), "patientwise", sep="_"), ".csv", sep="")
  }

  gcross_results = data.frame(matrix(ncol=radiusmax+2, nrow=0))

  XYdata$group = ""

  if (imagewise){
    XYdata$group = as.character(paste(XYdata[, patientcol], XYdata[, samplecol], sep=""))
  } else{
    XYdata$group = as.character(XYdata[,patientcol])
  }

  for (patient in unique(XYdata$group)){
    patient_data = XYdata[which(XYdata$group == patient),]
    cw <- convexhull.xy(patient_data[,xcol], patient_data[,ycol])
    ww <- owin(poly = cw$bdry)
    pp <- as.ppp(cbind(unlist(patient_data[,xcol]),unlist(patient_data[,ycol])), W = c(-1,2*max(patient_data[,xcol]),-1,2*max(patient_data[,ycol])))
    pp <- pp %mark% factor(unlist(patient_data$PastedPhenotype))
    pp$window <- ww

    if (sum(unlist(patient_data$PastedPhenotype) == "Reference") < 1 | sum(unlist(patient_data$PastedPhenotype) == "Non-Reference") < 1) {
      gcross_results = rbind.data.frame(gcross_results, c(patient, rep(0, radiusmax+1)), stringsAsFactors = FALSE)
    } else{
      Gcross1 = Gcross(pp,"Reference","Non-Reference",r = 0:radiusmax)
      gcross_results = rbind.data.frame(gcross_results, c(patient, t(Gcross1$km)), stringsAsFactors = FALSE)
    }
  }
  # write.table(gcross_results, "/Users/morganoneka/Documents/Grad School/Lab/Code/EAC_HGD/EAC/test.csv", sep=",", col.names=FALSE, quote=FALSE)

  write.table(gcross_results, gsub("//","/",fullfile(output, filename)), sep=",", col.names=FALSE, quote=FALSE)
}

colnames(summary_stats) <- c("CellType", "Count")
summary_stats <- distinct(summary_stats)
write.table(summary_stats, summaryfile, sep=",", col.names=TRUE, quote=FALSE, row.names=FALSE)
