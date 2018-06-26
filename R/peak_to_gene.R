
#each TSS is assigned its closest peak
#returns bed file with gene associated with closest peak to TSS
closestPeak = function(bed, outfile){

  file.remove(outfile)

  annot_bed = plyr::ddply(bed,plyr::.(chr),function(chr){
    annot = ChIPhandlr::ucsc[ChIPhandlr::ucsc$chrom == unique(chr$chr),]
    closest_idx = sapply(annot$TSS,function(x){
      return(which.min(abs(chr$peak_start-x)))
    })
    annot$bed_idx = closest_idx
    df = plyr::ddply(annot,plyr::.(bed_idx),function(sub){
      return(data.frame(bed_idx = unique(sub$bed_idx),
        hgnc_id = paste(unique(sub$hgnc_id),collapse = ",",sep = ""),
        ensembl_transcript_id = paste(unique(sub$ensembl_transcriptID), collapse = ",", sep = "")))
    })
    chr$hgnc_id = ""
    chr$ensembl_transcript_id = ""

    #fix factors issue
    chr[df$bed_idx,]$hgnc_id = as.character(df$hgnc_id)
    chr[df$bed_idx,]$ensembl_transcript_id = as.character(df$ensembl_transcript_id)

    write.table(chr,outfile,append = T,
      quote = F, col.names = F, row.names = F, sep = "\t")

    return("written")
  })

}


#each peak is assigned its closest TSS
#returns bed file containing all peaks in input bed annotated with with their closest TSS i.e. target_transcript, target_gene
closestTSS = function(bed, outfile){

  file.remove(outfile)

  annot_bed = plyr::ddply(bed,plyr::.(chr),function(chr){
    annot = ChIPhandlr::ucsc[ChIPhandlr::ucsc$chrom == unique(chr$chr),]
    closest_idx = sapply(chr$peak_start,function(x){
      return(which.min(abs(annot$TSS-x)))
    })

    annot = annot[closest_idx,]
    chr$target = annot$hgnc_id
    chr$target_transcript = annot$ensembl_transcriptID

    write.table(chr,outfile,append = T,
      quote = F, col.names = F, row.names = F, sep = "\t")


    return("written")
  })

}

#each TSS is assigned all peaks within some genomic window
#returns bed file annotated with the TSSs that the peak falls
#in the window for
fixedWindow = function(bed, outfile, window = c(-1000,1000)){

  file.remove(outfile)

  annot_ucsc = plyr::ddply(bed,plyr::.(chr),function(chr){
    annot = ChIPhandlr::ucsc[ChIPhandlr::ucsc$chrom == unique(chr$chr),]

    minus_strand = annot[annot$strand == "-",]
    plus_strand = annot[annot$strand == "+",]

    idx_plus = sapply(chr$peak_start,function(x){
      return(which(x > abs(plus_strand$TSS + window[1]) &
          x < abs(plus_strand$TSS + window[2])))
    })
    idx_minus = sapply(chr$peak_start,function(x){
      return(which(x < abs(minus_strand$TSS - window[1]) &
          x > abs(minus_strand$TSS - window[2])))
    })

    hgnc_ids =  paste(unlist(lapply(idx_plus,function(x){
      paste(unique(plus_strand[x,"hgnc_id"]),collapse = ",")})),
      unlist(lapply(idx_minus,function(x){
        paste(unique(minus_strand[x,"hgnc_id"]),collapse = ",")
      })),sep = ",")

    transcript_ids =  paste(unlist(lapply(idx_plus,function(x){
      paste(unique(plus_strand[x,"ensembl_transcriptID"]),collapse = ",")})),
      unlist(lapply(idx_minus,function(x){
        paste(unique(minus_strand[x,"ensembl_transcriptID"]),collapse = ",")
      })),sep = ",")

    hgnc_ids = gsub("^,*|,*$","",hgnc_ids)

    transcript_ids = gsub("^,*|,*$","",transcript_ids)

    chr$target_hgnc = hgnc_ids
    chr$target_tx = transcript_ids
    write.table(chr,outfile,append = T,
      quote = F, col.names = F, row.names = F, sep = "\t")
    return("written")
  })
}

#linear decay 1->0
linearPeakScore = function(bed, outfile, window = 50000){

  file.remove(outfile)

  annot = ChIPhandlr::ucsc
  annot = annot[annot$hgnc_id %in% genesetr::hgnc_dict_2018$approved,]
  minus_strand = annot[annot$strand == "-",]
  plus_strand = annot[annot$strand == "+",]

  #remove all but most upstream TSS for each gene
  minus_strand = minus_strand[order(minus_strand$TSS,decreasing = T),]
  plus_strand = plus_strand[order(plus_strand$TSS,decreasing = F),]

  minus_strand = minus_strand[!duplicated(minus_strand$hgnc_id),]
  plus_strand = plus_strand[!duplicated(plus_strand$hgnc_id),]
  annot = rbind(minus_strand, plus_strand)


  results = plyr::dlply(bed,plyr::.(TF),function(TF){

    chr_scores = return(plyr::ddply(TF,plyr::.(chr),function(chr){
      chr_annot = annot[annot$chrom == unique(chr$chr),]
      scores = sapply(chr_annot$TSS,function(tss){
        diff = abs(tss-chr$peak_start)
        diff = diff[diff<window]
        scores = 1-diff/window
        return(sum(scores))})
      return(data.frame(gene = chr_annot$hgnc_id, score = scores))
    }))

    chr_scores$TF = unique(TF$TF)
    write.table(chr_scores,outfile,append = T,
      quote = F, col.names = F, row.names = F, sep = "\t")
    return("written")

  })
}




#exponential decay 1->0
exponentialPeakScore = function(bed, outfile, window = c(-100000,100000)){

}


