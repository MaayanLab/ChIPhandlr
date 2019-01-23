
#' @export

closestTSS = function(bed, chr_col, peak_col, outfile){

  file.remove(outfile)
  annot = ChIPhandlr::ucsc
  annot = annot[annot$hgnc_id %in% unique(genesetr::hgnc_dict),]

  annot_bed = plyr::ddply(bed,chr_col,function(chr){
    annot_chr = annot[annot$chrom == unique(chr[,chr_col]),]
    closest_idx = sapply(chr[,peak_col], function(x){
      return(which.min(abs(annot$TSS-x)))
    })

    annot_close = annot_chr[closest_idx, ]
    chr$target = annot_close$hgnc_id
    chr$target_transcript = annot_close$ensembl_transcriptID
    chr$TSS = annot_close$TSS

    write.table(chr,outfile,append = T,
      quote = F, col.names = F, row.names = F, sep = "\t")


    return("written")
  })

}
