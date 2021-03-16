library("optparse")
library("RJSONIO")

option_list = list( make_option(c("-i", "--input"), type="character", default=NULL,
                                help="Input file or directory"),
                    make_option(c("-o", "--output"), type="character", default=NULL,
                                help="Output file"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
input = opt$input
output = opt$output

phenotypes <- fromJSON(input)

filenames = c()
for (inx in phenotypes){
  ref = inx[[1]][[1]]$Reference
  ref_phenotypes = ""
  for (option in ref){
    ref_phenotypes = paste(ref_phenotypes, paste(names(option), option, sep="_", collapse="_"), sep="_or_")
  }
  ref_phenotypes = str_split_fixed(ref_phenotypes, "_", 3)[3]
  
  nonref = inx[[1]][[1]]$`Non-Reference`
  nonref_phenotypes = ""
  for (option in nonref){
    # paste(names(option), option, sep="_", collapse="_")
    nonref_phenotypes = paste(nonref_phenotypes, paste(names(option), option, sep="_", collapse="_"), sep="_or_")
  }
  nonref_phenotypes = str_split_fixed(nonref_phenotypes, "_", 3)[3]
  
  filenames = c(filenames, paste(ref_phenotypes, nonref_phenotypes, sep="_VS_"))
}

# x <- unlist(phenotypes)
# paste(str_split_fixed(names(x), "\\.",3)[,3], x, sep="_")

write.table(filenames, output, col.names=FALSE, row.names=FALSE, quote=FALSE)
