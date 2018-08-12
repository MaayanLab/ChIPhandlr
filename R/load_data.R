#returns MACS xls object
loadMACSxls = function(filename){
  df = read.table(filename, sep = "\t", header = T,
    comment.char = "#", stringsAsFactors = F, quote="")
  colnames(df) = as.character(df[1,])
  df = df[-1,]
  return(df)
}