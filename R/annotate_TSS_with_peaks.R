#takes in a list of beds from multiple experiments
#returns annotated TSSs where each TSS is annotated with its closest peak from among the bed list
#can accomodate one bed or a list of beds in data frame format

#' @useDynLib ChIPhandlr
#' @importFrom Rcpp sourceCpp
NULL
#' @export

closestPeak = function(bed_list, peak_col, chr_col = "chr", exp_col, dir = "both"){
  #stack bed data frames
  bed = plyr::rbind.fill(bed_list)
  annot = ChIPhandlr::ucsc
  annot = annot[annot$hgnc_id %in% unique(genesetr::hgnc_dict),]

  #slice bed file by chromosome
  annot_TSS = plyr::ddply(bed,chr_col,function(chr){


    #subset TSS data frame by chromosome
    annot_chr = annot[annot$chrom == unique(chr[,chr_col]),]

    #determine bed index of closest peak
    if(dir == "both"){
      closest_idx = sapply(annot_chr$TSS,function(x){
        return(which.min(abs(chr[,peak_col]-x)))})

      #annotate TSSs with closest peak
      annot_chr$peak_name = as.character(chr[closest_idx,exp_col])
      annot_chr$peak_loc = chr[closest_idx,peak_col]

    }else{
      plus_strand = annot_chr[annot_chr$strand == "+",]
      minus_strand = annot_chr[annot_chr$strand == "-",]
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



        #annotate TSSs with closest peak
        annot_chr$peak_name = NA
        annot_chr$peak_loc = NA

        if(length(closest_idx_plus) == 0) closest_idx_plus = NA
        if(length(closest_idx_minus) == 0) closest_idx_minus = NA

        annot_chr[annot_chr$strand == "+",
          "peak_name"] = as.character(chr[closest_idx_plus, exp_col])
        annot_chr[annot_chr$strand == "+",
          "peak_loc"] = chr[closest_idx_plus, peak_col]

        annot_chr[annot_chr$strand == "-",
          "peak_name"] = as.character(chr[closest_idx_minus, exp_col])
        annot_chr[annot_chr$strand == "-",
          "peak_loc"] = chr[closest_idx_minus, peak_col]

      }else{
        error("dir argument must be 'upstream', 'downstream' or 'both'")
      }

    }

    return(annot_chr)
    })

  return(annot_TSS)

} #end function

#each TSS is assigned all peaks within some genomic window
#returns TSSs annotated with all peaks that fall within the genomic window
#' @export

fixedWindow = function(bed_list, outfile, peak_col, chr_col = "chr", exp_col, window = c(-1000,1000)){

  file.remove(outfile)
  bed = plyr::rbind.fill(bed_list)

  annot = ChIPhandlr::ucsc
  annot = annot[annot$hgnc_id %in% unique(genesetr::hgnc_dict),]

  annot_TSS = plyr::ddply(bed, chr_col, function(chr){

    annot_chr = annot[annot$chrom == unique(chr[,..chr_col]), ]

    #separate UCSC annotations into plus and minus strand data frames
    minus_strand = annot_chr[annot_chr$strand == "-",]
    plus_strand = annot_chr[annot_chr$strand == "+",]

    #generate a list with one entry per TSS. Each entry of the list contains the indices of
    #the peaks that fall within the window of that TSS.
    idx_plus = lapply(plus_strand[,"TSS"], function(tss){
      return(which(window[1] <= (chr[,..peak_col] - tss) &
          window[2] >= (chr[,..peak_col] - tss)))
    })

    #generate a list with one entry per TSS. Each entry of the list contains the indices of
    #the peaks that fall within the window of that TSS.
    idx_minus = lapply(minus_strand[1:100,"TSS"], function(tss){
      return(which(window[1] <= (tss - chr[,..peak_col]) &
          window[2] >= (tss - chr[,..peak_col])))
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
#' @export

peakScoreDist = function(bed, score_col, peak_col, chr_col = "chr", window = 20000, decay = "linear", d = NULL){
  ptm = proc.time()


  #auxiliary function to restrict base diag
  diag_mat = function(x){
    if(length(x)>1){
      return(diag(x))
    }else{
      return(x)
    }
  }

  #internal package annotations (at this point, only hg38 available)
  annot = ChIPhandlr::ucsc
  annot = annot[annot$hgnc_id %in% unique(genesetr::hgnc_dict),]


  #slice bed into chromosomes
  results = plyr::ddply(bed,plyr::.(chr),function(chr){


    #subset ucsc annotations for chromosome
    chr_annot = annot[annot$chrom == unique(chr[,chr_col]),]


    results_sub = data.frame(ensembl = chr_annot$ensembl_transcriptID,
      tss = chr_annot$TSS,
      strand = chr_annot$strand,
      gene = chr_annot$hgnc_id,
      stringsAsFactors = F)


    #compute distance matrix between all peak locations on
    #on chromosome and all TSSs on chromosome
    diff = abs(outer(chr_annot$TSS,as.numeric(chr[,peak_col]),'-'))
    peak_scores = as.numeric(chr[,score_col])

    #scoring function decays linearly with distance from TSS
    if(decay == "linear"){
      #multiple bed peak scores by distance to TSS where scaling
      #coefficient is 1 when there is 0 distance between the peak
      #and the TSS and 0 when the distance between the peak and
      #the TSS reaches the edge of the TSS w
      ptm = proc.time()
      for(i in 1:length(window)){
        scored_diff = 1-diff/window[i]
        for(j in 1:length(peak_scores)){
          scored_diff[,j] = scored_diff[,j]*peak_scores[j]
        }

        #negative scores mean the peak is outside of the desired
        #TSS genomic window; set these scores to 0
        scored_diff = replace(scored_diff, scored_diff<0, 0)
        results_sub[,paste("window_",window[i]/1000,"k",sep = "")] = rowSums(scored_diff)
      }




      #scoring function decays exponentially with distance from TSS
    }else if(decay == "exp"){
      if(length(window)>1) error("May only choose one window when using exponential decay. May have multiple decay constants.")

      #for each decay constant
      for(i in 1:length(d)){
        scored_diff = exp(-diff/d[i])
        for(j in 1:length(peak_scores)){
          scored_diff[,j] = scored_diff[,j]*peak_scores[j]
        }

        scored_diff = replace(scored_diff, scored_diff<0, 0)
        results_sub[,paste("decayconst_",d[i],sep = "")] = rowSums(scored_diff)
      }
    }else{
      error("Argument decay must have values 'linear' or 'exp'.")
    }

    #return data frame with all TSSs and associated score
    #corresponding to sum of scaled
    #peak intensities within designated genomic window around
    #each TSS
    return(results_sub)
  }) #end aggregation of dataframes for all chromosomes
      #into one data frame containing
      #all annotated TSSs
  proc.time()-ptm

  return(results)

}#end function peakScoreDist()

