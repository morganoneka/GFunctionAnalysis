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

getPhenotypeString <- function(ref, verbose){
  ref_phenotypes = ""
  for (option in ref){
    if (verbose){
      ref_phenotypes = paste(ref_phenotypes, paste(names(option), option, sep="_", collapse="_"), sep="_or_")  
    } else{
      ref_phenotypes = paste(ref_phenotypes, paste(option, sep="_", collapse="_"), sep="_or_")
    }
    
  }
  return (str_split_fixed(ref_phenotypes, "_", 3)[3])
}

getRowsWithPhenotype <- function(ref, XYdata, verbose){
  ref_indices = c()
  for (def in ref){
    ref_cols = sort(names(def))
    
    phenotype_string = ""
    
    if (verbose){
      phenotype_string = paste(rbind(ref_cols, sort(unlist(def))), collapse="_")  
    } else{
      phenotype_string = paste(sort(unlist(def)), sep="_", collapse="_")
    }
    
    if (verbose){
      tmp <- zipFastener(data.frame( matrix(sort(rep(ref_cols, nrow(XYdata))), ncol=length(ref_cols)) ), data.frame(XYdata[,ref_cols]))
      XYdata$ref <- apply( tmp, 1 , paste , collapse = "_" )
    } else{
      if (length(ref_cols) > 1){
        XYdata$ref = apply( XYdata[,sort(unlist(ref_cols))], 1, paste, collapse="_")  
      } else{
        XYdata$ref = XYdata[,ref_cols]
      } 
      
    }
    
    ref_indices = c(ref_indices, which(XYdata$ref == phenotype_string))
  }
  
  return(unique(ref_indices))
}
