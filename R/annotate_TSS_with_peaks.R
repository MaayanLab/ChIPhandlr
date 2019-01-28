#takes in a list of beds from multiple experiments
#returns annotated TSSs where each TSS is annotated with its closest peak from among the bed list
#can accomodate one bed or a list of beds in data frame format


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
#' @import magrittr

fixedWindow = function(bed_list, peak_col, chr_col = "chr", exp_col, window = c(-1000,1000)){


  bed = plyr::rbind.fill(bed_list)

  annot = ChIPhandlr::ucsc
  annot = annot[annot$hgnc_id %in% unique(genesetr::hgnc_dict),]

  fixed_win = plyr::dlply(bed, chr_col, function(chr){

    annot_chr = annot[annot$chrom == unique(unlist(chr[,chr_col])), ]

    #separate UCSC annotations into plus and minus strand data frames
    minus_strand = annot_chr[annot_chr$strand == "-",]
    plus_strand = annot_chr[annot_chr$strand == "+",]

    #take only the first TSS for each gene symbol
    minus_strand = minus_strand[order(-minus_strand$TSS),]
    minus_strand = minus_strand[!duplicated(minus_strand$hgnc_id),]
    plus_strand = plus_strand[order(plus_strand$TSS),]
    plus_strand = plus_strand[!duplicated(plus_strand$hgnc_id),]


    #separate by TF
    tf_results = plyr::dlply(chr, "TF", function(tf){
      plus_diff = outer(plus_strand$TSS,as.numeric(tf$peak_start),'-')
      plus_diff_log = -plus_diff > window[1] & -plus_diff < window[2]
      plus_targets = unique(plus_strand[rowSums(plus_diff_log) > 0,"hgnc_id"])

      rm(plus_diff, plus_diff_log)

      minus_diff = abs(outer(minus_strand$TSS,as.numeric(tf$peak_start),'-'))
      minus_diff_log = minus_diff > window[1] & minus_diff < window[2]
      minus_targets = unique(minus_strand[rowSums(minus_diff_log) > 0,"hgnc_id"])

      return(unique(c(minus_targets,plus_targets)))

    })

   return(tf_results)
  })

  #combine tf_results
  all_tfs = unique(bed$TF)

  results_gmt = lapply(all_tfs,function(tf){
    return(lapply(fixed_win,function(chr_results){
      if(tf %in% names(chr_results)){
        return(chr_results[[tf]])
      }
    }) %>% unlist %>% as.character %>% unique)
  })
  names(results_gmt) = all_tfs
  return(results_gmt)

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

    #further slice bed into TFs


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

