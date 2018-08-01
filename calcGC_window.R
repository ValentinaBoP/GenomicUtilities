calcGC_window = function(file, fasta = NULL, window, pattern = character()){
  
  # Valentina Peona 01/07/2018
  
  # usage:
  # GCcontent_df = calcGC_window(file = "genome.fasta", window = 10000)
  
  # required libraries
  require(Biostrings)
  require(BSgenome)
  
  # check window size
  if(length(window) == 0 | window <= 0){
    
    window = 10000
    
  }
  
  # check if the correct data have been provided
  if(length(fasta) == 0){
    
    print("reading fasta file")
    fasta = readDNAStringSet(filepath = file)
    print("fasta file imported")
    
  }
  
  # if pattern has been given, subset the sequences and names by it
  if(length(pattern) > 0){
    
    boo = which(fasta@ranges@NAMES %in% pattern)
    fasta = fasta[boo]
    
  }
  
  chrWidth = data.frame(chrom = fasta@ranges@NAMES, width = fasta@ranges@width)
  
  GCcontent = createWindows(data = chrWidth, window = window)
  
  print("calculating GC content")
  list_seqs = strsplit(as.character((getSeq(fasta, as(GCcontent, "GRanges")))), split = "")
  GCcontent$gc = sapply(list_seqs, GC)
  
  return(GCcontent)
  
}

createWindows = function(data, window){
  
  newData = data.frame(chrom = character(), start = integer(), end = integer())
  
  for(i in 1:nrow(data)){
    
    start = as.integer(seq(from = 1, to = data$width[i], by = window))
    end = as.integer(c(start[2:length(start)] - 1, data$width[i]))
    
    newData = rbind(newData, data.frame(chrom = data$chrom[i], start = start, end = end))
    
  }
  
  return(newData)
  
}
