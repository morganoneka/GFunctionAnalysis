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
  make_option(c("-r", "--radius"), type="character", default="60,120,200",
              help="Radii to calculate AUC for")
)

# parse the command line arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input = opt$input
output = opt$output
radius = opt$radius

# get a list of all files in the input directory
files <- list.files(path=input, pattern="*_KM.csv", full.names=FALSE, recursive=FALSE)

all_radii <- as.numeric(strsplit(radius, ",")[[1]])

auc_data <- data.frame(matrix(nrow=0, ncol= (length(files)*length(all_radii) + 1)))
colnames(auc_data) <- c("sample", (1:(ncol(auc_data)-1)))

col_strings <- c("sample")

for (inputfile in files){
  # get interaction from filename
  interaction <- strsplit(inputfile, "\\.csv")[[1]][1]
  
  for (i in length(all_radii)){
    col_strings <- c(col_strings, paste(interaction, all_radii, sep="_"))
  }
  
  # read in file
  GcrossData <- read.table(fullfile(input,inputfile), header = FALSE, sep = ",", fill = TRUE)
  
  # create sub_df for just this interaction
  interaction_auc <- data.frame(matrix(nrow=0, ncol = (1+all_radii)))
  
  # each row = 1 patient/sample
  for (row in 1:nrow(GcrossData)){
    line <- GcrossData[row,]
    
    # first column = row, second column = sample
    sample_no <- line[,2]
    
    # get data points
    gcross_points <- line[,3:length(line)]
    
    single_row_data <- c(sample_no)
    
    # calculate auc 
    for (rad in all_radii){
      # calculate auc up to rad+1 because the first data point is 0
      auc <- trapz(0:rad, t(gcross_points[1:(rad+1)]))
      single_row_data <- cbind(single_row_data, auc)
    }
    
    # add this sample's data to dataframe for all interactions
    interaction_auc <- rbind(interaction_auc, single_row_data)
  }
  
  colnames(interaction_auc)[1] <- ("sample")
  
  if (nrow(auc_data) == 0){
    auc_data = interaction_auc
  } else{
    auc_data <- merge(auc_data, interaction_auc, by="sample")  
  }
  
  colnames(auc_data) <- col_strings
  
}

write.table(auc_data,fullfile(output, "GcrossAUC.csv"),  sep=",", col.names=TRUE, row.names=FALSE)
