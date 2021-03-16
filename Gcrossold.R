library("spatstat")
library("matlab")
library("optparse")
library("RJSONIO")
library("dplyr")

# zipFastener for TWO dataframes of unequal length
zipFastener <- function(df1, df2, along=2)
{
  # parameter checking
  if(!is.element(along, c(1,2))){
    stop("along must be 1 or 2 for rows and columns
         respectively")
  }
  # if merged by using zip feeding along the columns, the
  # same no. of rows is required and vice versa
  if(along==1 & (ncol(df1)!= ncol(df2))) {
    stop ("the no. of columns has to be equal to merge
          them by zip feeding")
  }
  if(along==2 & (nrow(df1)!= nrow(df2))) {
    stop ("the no. of rows has to be equal to merge them by
          zip feeding")
  }
  
  # zip fastener preperations
  d1 <- dim(df1)[along]
  d2 <- dim(df2)[along]
  i1 <- 1:d1           # index vector 1
  i2 <- 1:d2 + d1      # index vector 2
  
  # set biggest dimension dMax
  if(d1==d2) {
    dMax <- d1
  } else if (d1 > d2) {
    length(i2) <- length(i1)    # make vectors same length, 
    dMax <- d1                  # fill blanks with NAs   
  } else  if(d1 < d2){
    length(i1) <- length(i2)    # make vectors same length,
    dMax <- d2                  # fill blanks with NAs   
  }
  
  # zip fastener operations
  index <- as.vector(matrix(c(i1, i2), ncol=dMax, byrow=T))
  index <- index[!is.na(index)]         # remove NAs
  
  if(along==1){
    colnames(df2) <- colnames(df1)   # keep 1st colnames                  
    res <- rbind(df1,df2)[ index, ]  # reorder data frame
  }
  if(along==2) res <- cbind(df1,df2)[ , index]           
  
  return(res)
  }

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
              help="Output file"),
  make_option(c("-j", "--jsonphenotype"), type="character", default=NULL,
              help="Json file containing cell phenotypes"),
  make_option(c("--sep"), type="character", default=",",
              help="Separator used in input data")
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

# read in XYdata, phenotype info
XYdata <- read.table(input, header = TRUE, sep = sep, fill = TRUE)
print(head(XYdata))
phenotypes <- fromJSON(jsonphenotype)

# TODO: delete me, i'm just here for testing
# phenotypes <- list()
# phenotypes[[1]] <- list(Epithelial = "pos")
# phenotypes[[2]] <- list(CD8 = "pos", Tcell = "pos")

# identify reference and non-reference cells
ref_cols = names(unlist(phenotypes[[1]]))
ref_vals = unlist(phenotypes[[1]])
ref_phenotype = paste(rbind(ref_cols, ref_vals), collapse="_")

nonref_cols = names(unlist(phenotypes[[2]]))
nonref_vals = unlist(phenotypes[[2]])
nonref_phenotype = paste(rbind(nonref_cols, nonref_vals), collapse="_")

tmp <- zipFastener(data.frame( matrix(sort(rep(nonref_cols, nrow(XYdata))), ncol=length(nonref_cols)) ), data.frame(XYdata[,sort(nonref_cols)]))
XYdata$nonref <- apply( tmp, 1 , paste , collapse = "_" )
tmp <- zipFastener(data.frame( matrix(sort(rep(ref_cols, nrow(XYdata))), ncol=length(ref_cols)) ), data.frame(XYdata[,sort(ref_cols)]))
XYdata$ref <- apply( tmp, 1 , paste , collapse = "_" )

XYdata$PastedPhenotype = ""
XYdata[which(XYdata$nonref == nonref_phenotype), "PastedPhenotype"] = "Non-Reference"
XYdata[which(XYdata$ref == ref_phenotype), "PastedPhenotype"] = "Reference"

# create groups based on patient, or patient and sample
gcross_results = data.frame(matrix(ncol=radiusmax+2, nrow=0))
for (patient in unique(XYdata[,patientcol])){
  patient_data = XYdata[which(XYdata[,patientcol] == patient),]
  cw <- convexhull.xy(patient_data[,xcol], patient_data[,ycol])
  ww <- owin(poly = cw$bdry)
  pp <- as.ppp(cbind(unlist(patient_data[,xcol]),unlist(patient_data[,ycol])), W = c(-1,2*max(patient_data[,xcol]),-1,2*max(patient_data[,ycol])))
  pp <- pp %mark% factor(unlist(patient_data$PastedPhenotype))
  pp$window <- ww
  
  if (sum(unlist(patient_data$PastedPhenotype) == "Reference") < 1 | sum(unlist(patient_data$PastedPhenotype) == "Non-Reference") < 1) {
    gcross_results = rbind(gcross_results, c(patient, rep(0, radiusmax+1)))
  } else{
    Gcross1 = Gcross(pp,"Reference","Non-Reference",r = 0:radiusmax)
    gcross_results = rbind(gcross_results, c(patient, t(Gcross1$km)))
  }
}

write.table(gcross_results, output, sep=",", col.names=FALSE)
