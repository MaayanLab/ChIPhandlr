#returns MACS xls object
loadMACSxls = function(filename){
  header_line = 28
  return(read.table(filename, sep = "\t", header = T, skip = header_line-1))
}