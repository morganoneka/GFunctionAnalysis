library("spatstat")
library("matlab")
library("optparse")
library("ggplot2")
library("pracma")
library("stringr")
source("/Users/morganoneka/Box/My Stuff/Gcross/functions.R")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input file or directory"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory"),
  make_option(c("-s", "--step"), type="numeric", default=10,
              help="Step size to calculate AUC for"),
  make_option(c("-e", "--endpoint"), type="numeric", default=200,
              help="The end point to stop calculating bins")
)


# parse the command line arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input = opt$input
output = opt$output
step = as.numeric(opt$step)
endpoint = as.numeric(opt$endpoint)

# get a list of all files in the input directory
files <- list.files(path=input, pattern="*_KM.csv", full.names=FALSE, recursive=FALSE)

# bin start points
sequences <- seq(0, endpoint, step)

for (inputfile in files){
  auc_data <- data.frame(matrix(nrow=0, ncol= (length(sequence) + 1)))
  colnames(auc_data) <- c("sample", (1:(ncol(auc_data)-1)))
  
  col_strings <- c("sample")
  for (i in 1:(length(sequences)-1)){
    bin_start = sequences[i]
    bin_end = sequences[i+1]
    
    col_strings <- c(col_strings, paste(bin_start,bin_end,sep="_"))
  }
  
  print(col_strings)
  
  # get interaction from filename
  interaction <- strsplit(inputfile, "\\.csv")[[1]][1]
  
  # read in data
  GcrossData <- read.table(fullfile(input,inputfile), header = FALSE, sep = ",", fill = TRUE)
  
  # iterate over all individuals
  for (row in 1:nrow(GcrossData)){
    line <- GcrossData[row,]
    
    # first column = row, second column = sample
    sample_no <- line[,2]
    
    # get data points
    gcross_points <- line[,3:length(line)]
    
    cumArea = cumtrapz(sequences,unlist(gcross_points[sequences + 2]))
    vec_AUCdiffs = cumArea[2:length(sequences)] - cumArea[1:(length(sequences)-1)]
    
    auc_data <- rbind(auc_data, c(sample_no, vec_AUCdiffs), stringsAsFactors=FALSE)
  }
  
  output_filename = paste("AUC_Kbins_", interaction, ".csv", sep="")
  colnames(auc_data) <- col_strings
  write.table(auc_data,fullfile(output, output_filename),  sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
}