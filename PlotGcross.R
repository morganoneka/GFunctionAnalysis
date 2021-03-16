library("optparse")
library("ggplot2")
library("matlab")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input file or directory"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory"),
  make_option(c("--summary"), type="character", default=NULL,
              help="Name for file to save gcross summary info to"),
  make_option(c("-w", "--imagewise"), default=FALSE, help="Use this flag for image-wise gcross")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
input = opt$input
output = opt$output
summaryfile = opt$summary
imagewise = opt$imagewise

files <- list.files(path=input, pattern="*_KM.csv", full.names=FALSE, recursive=FALSE)

for (inputfile in files){
  XYdata <- read.table(fullfile(input,inputfile), header = FALSE, sep = ",", fill = TRUE, stringsAsFactors=FALSE)
  print(inputfile)
  print(nrow(XYdata))
  for (i in 1:nrow(XYdata)){
    row = XYdata[i,3:ncol(XYdata)]
    coordinates = data.frame(x=0:(length(row)-1), y=t(row))
    colnames(coordinates) <- c("x", "y")
    ggplot(coordinates, aes(x=x, y=y)) + geom_line(size=1) + ggtitle(paste("G-cross curve for", XYdata[i,2]) )
    filename=""
    sample<-strsplit(tail(strsplit(XYdata[i,2], "/")[[1]],1), "\\.")[[1]][1]
    if (imagewise){
      filename = paste(paste("sample", sample, strsplit(inputfile, "\\.")[[1]][1], "imagewise", sep="_"), ".png", sep="")
    } else{
      filename = paste(paste("sample", sample, strsplit(inputfile, "\\.")[[1]][1], "patientwise", sep="_"), ".png", sep="")
    }
    print(filename)
    ggsave(fullfile(output,filename), plot = last_plot(), device="png")
  }
}



summary_stats = data.frame(matrix(ncol=0, nrow=0))
write.table(summary_stats, summaryfile, sep=",", col.names=TRUE)