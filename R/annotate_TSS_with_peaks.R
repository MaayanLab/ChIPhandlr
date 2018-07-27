#takes in a list of beds from multiple experiments
#returns annotated TSSs where each TSS is annotated with its closest peak from among the bed list
#can accomodate one bed or a list of beds in data frame format
closestPeak = function(bed_list, peak_col, chr_col = "chr", exp_col, dir = "both"){
  #stack bed data frames
  bed = plyr::rbind.fill(bed_list)

  #slice bed file by chromosome
  annot_TSS = plyr::ddply(bed,chr_col,function(chr){
    print(unique(chr$chr))

    #subset TSS data frame by chromosome
    annot = ChIPhandlr::ucsc[ChIPhandlr::ucsc$chrom == unique(chr[,chr_col]),]

    #determine bed index of closest peak
    if(dir == "both"){

      closest_idx = sapply(annot$TSS,function(x){
        return(which.min(abs(chr[,peak_col]-x)))})

      #annotate TSSs with closest peak
      annot$peak_name = as.character(chr[closest_idx,exp_col])
      annot$peak_loc = chr[closest_idx,peak_col]

    }else{
      plus_strand = annot[annot$strand == "+",]
      minus_strand = annot[annot$strand == "-",]

      if(dir == "upstream"){
        closest_idx_plus = sapply(plus_strand$TSS,function(x){
          diff = x - chr[,peak_col]
          diff[diff<0] = NA
          if(any(!is.na(diff))){
            return(which.min(diff))
          }else return(NA)})

        closest_idx_minus = sapply(minus_strand$TSS,function(x){
          diff = chr[,peak_col] - x
          diff[diff<0] = Inf
          if(any(!is.na(diff))){
            return(which.min(diff))
          }else return(NA)})

      }else if(dir == "downstream"){
        closest_idx_plus = sapply(plus_strand$TSS,function(x){
          diff = chr[,peak_col] - x
          diff[diff<0] = NA
          if(any(!is.na(diff))){
            return(which.min(diff))
          }else return(NA)})

        closest_idx_minus = sapply(minus_strand$TSS,function(x){
          diff = x - chr[,peak_col]
          diff[diff<0] = Inf
          if(any(!is.na(diff))){
            return(which.min(diff))
          }else return(NA)})
      }else{
        error("dir argument must be 'upstream', 'downstream' or 'both'")
      }

      #annotate TSSs with closest peak
      annot$peak_name = NA
      annot$peak_loc = NA

      if(length(closest_idx_plus) == 0) closest_idx_plus = NA
      if(length(closest_idx_minus) == 0) closest_idx_minus = NA

      annot[annot$strand == "+", "peak_name"] = as.character(chr[closest_idx_plus, exp_col])
      annot[annot$strand == "+", "peak_loc"] = chr[closest_idx_plus, peak_col]

      annot[annot$strand == "-", "peak_name"] = as.character(chr[closest_idx_minus, exp_col])
      annot[annot$strand == "-", "peak_loc"] = chr[closest_idx_minus, peak_col]

    }

    return(annot)})

  return(annot_TSS)

} #end function

#each TSS is assigned all peaks within some genomic window
#returns TSSs annotated with all peaks that fall within the genomic window
#can accomodate one bed or a list of beds in data frame format
fixedWindow = function(bed_list, peak_col, chr_col = "chr", exp_col, window = c(-1000,1000)){


  annot_TSS = plyr::ddply(bed, chr_col, function(chr){
    annot = ChIPhandlr::ucsc[ChIPhandlr::ucsc$chrom == unique(chr$chr), ]

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
#returns annotated TSS file
#annotates all transcription start sites with score
#can accomodate one bed in data frame format

linearPeakScore = function(bed_list, score_col, peak_col, chr_col = "chr", window = 20000, decay = "linear", b = NULL){

  annot = ChIPhandlr::ucsc

  scored_TSSs = plyr::ddply(bed,plyr::.(chr),function(chr){

    chr_annot = annot[annot$chrom == unique(chr[,chr_col]),]
    diff = abs(outer(chr_annot$TSS,chr[,peak_col],'-'))

    if(decay == "linear"){
      scores = (1-diff/window) %*% diag(chr[,score_col])
      scores[scores<0] = 0
    }else if(decay == "exp"){
      scores = (exp(-diff)) %*% diag(chr[,score_col])

    }else{
      error("Argument decay must have values 'linear' or 'exp'.")
    }

    return(data.frame(ensembl = chr_annot$ensembl_transcriptID,
      tss = chr_annot$TSS,
      strand = chr_annot$strand,
      gene = chr_annot$hgnc_id,
      score = rowSums(scores),
      stringsAsFactors = F))
  })
  return(scored_TSSs)
}
