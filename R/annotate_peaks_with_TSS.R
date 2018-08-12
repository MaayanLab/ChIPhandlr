closestTSS = function(bed, outfile){

  file.remove(outfile)

  annot_bed = plyr::ddply(bed,plyr::.(chr),function(chr){
    annot = ChIPhandlr::ucsc[ChIPhandlr::ucsc$chrom == unique(chr$chr),]
    closest_idx = sapply(chr$peak_start, function(x){
      return(which.min(abs(annot$TSS-x)))
    })

    annot = annot[closest_idx, ]
    chr$target = annot$hgnc_id
    chr$target_transcript = annot$ensembl_transcriptID

    write.table(chr,outfile,append = T,
      quote = F, col.names = F, row.names = F, sep = "\t")


    return("written")
  })

}
