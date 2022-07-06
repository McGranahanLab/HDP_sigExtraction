#######################################################
########      YAPSA indel_functions.R script   ########
#######################################################
#copied from: https://rdrr.io/bioc/YAPSA/src/R/indel_functions.R on the 30.09.2021 at 17:30

# README.md
# 1. Usage of YAPSA
# 2. Signature-specific cutoffs
# 3. Confidence Intervals
# 4. Stratified Analysis of Mutational Signatures
# 5. Indel signature analysis
# 6. Usage of YAPSA for Whole Exome Sequencing (WES) data
# Home / Bioconductor / YAPSA / R/indel_functions.R
# ads via Carbon
# Code in Node.js, Java, Python and other open-source languages.
# ADS VIA CARBON
# R/indel_functions.R
# In YAPSA: Yet Another Package for Signature Analysis
# Defines functions plotExposuresConfidence_indel confidence_indel_only_calulation confidence_indel_calulation create_indel_mutation_catalogue_from_df create_indel_mut_cat_from_df attribution_of_indels attribute_sequence_contex_indel getSequenceContext
# Documented in attribute_sequence_contex_indel attribution_of_indels confidence_indel_calulation confidence_indel_only_calulation create_indel_mutation_catalogue_from_df create_indel_mut_cat_from_df getSequenceContext plotExposuresConfidence_indel
# Copyright © 2014-2019  The YAPSA package contributors
# This file is part of the YAPSA package. The YAPSA package is licenced under
# GPL-3

#' Extracts the sequence context up and downstream of a nucleotide position
#'
#' @param position Start position of the considered INDEL
#' @param chr Chromosome of the considered INDEL
#' @param offsetL Number of nucleotides downstream of \code{position}
#' @param offsetR Number of nucleotides upstream of \code{position}
#'
#' @return Returns a character string containing the defined seqeunce context
#'
#' @export
#' @importFrom Biostrings getSeq DNAString subseq
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @examples
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' sequence_context <- getSequenceContext(position = 123456789, chr = "chr12",
#'                                        offsetL= 10, offsetR=50)
#' sequence_context
#' 
getSequenceContext <- function(position, chr, offsetL= 10, offsetR=50){
  if(is.numeric(chr)){
    chr <- paste0("chr", chr)
    sequence <- getSeq(Hsapiens, chr, position-offsetL, position+offsetR)
    context_string <- as.character(sequence)
    sequenceContext <- list(sequence=sequence,context_string=context_string)
    return(sequenceContext)
  }else{
    sequence <- getSeq(Hsapiens, chr, position-offsetL, position+offsetR)
    context_string <- as.character(sequence)
    sequenceContext <- list(sequence=sequence,context_string=context_string)
    return(sequenceContext)
  }
}


#' Attribution of sequence context and size for an INDEL
#'
#' The function is a wrapper and uses \code{\link[YAPSA]{getSequenceContext}}
#' to annotate the sequence context.
#'
#' @param in_dat VRanges object or data frame which carries one column for the
#'   reference base and one column for the variant base
#' @param in_REF.field String indicating which column of \code{in_dat} carries
#'   the reference base if dealing with data frames
#' @param in_ALT.field String indicating which column of \code{in_dat} carries
#'   the variant base if dealing with data frames
#' @param in_verbose Verbose if \code{in_verbose=1}
#' @param in_offsetL Number of nucleotides which should be annotated downstream
#'   of the variant. Per default 10 bps are annotated
#' @param in_offsetR Number of nucleotides which should be annotated upstream of
#'   the catiant. Per default 50 bps are annotated
#'
#' @return VRanges object or data frame with the same number rows and additional
#' columns containing the type of INDEL (Ins = insertion and Del = deletion),
#' the annotated sequence context of the defined length, the absolute number of
#' exchanged nucleotides and the nucleotide exchange between \code{in_REF.field}
#' and \code{in_ALT.field}.
#'
#' @export
#'
#' @examples
#' data(GenomeOfNl_raw)
#' GenomeOfNl_context <- attribute_sequence_contex_indel(
#'                                    in_dat = head(GenomeOfNl_raw),
#'                                    in_REF.field = "REF",
#'                                    in_ALT.field = "ALT",
#'                                    in_verbose = FALSE,
#'                                    in_offsetL= 10, in_offsetR=50)
#' GenomeOfNl_context
#' 
attribute_sequence_contex_indel <- function(in_dat, in_REF.field = "REF",
                                            in_ALT.field = "ALT",
                                            in_verbose = FALSE, 
                                            in_offsetL = 10,
                                            in_offsetR = 50) {
  print(paste("Indel sequence context attribution of total ", dim(in_dat)[1],
              " indels. This could take a while..."))
  name_list <- names(in_dat)
  ## exception handling for input fields
  if(tolower(in_REF.field) %in% tolower(name_list)){
    if(in_verbose) cat("YAPSA:::attribute_nucleotide_exchanges::in_REF.field",
                       " found. Retrieving REF information.\n")
    REF_ind <- min(which(tolower(name_list)==tolower(in_REF.field)))
  }else{
    if(in_verbose) cat("YAPSA:::attribute_nucleotide_exchanges::error: ",
                       "in_REF.field not found. Return NULL.\n")
    return(NULL)
  }
  if(tolower(in_ALT.field) %in% tolower(name_list)) {
    if(in_verbose) cat("YAPSA:::attribute_nucleotide_exchanges::in_ALT.field",
                       " found. Retrieving ALT information.\n")
    ALT_ind <- min(which(tolower(name_list)==tolower(in_ALT.field)))
  }else{
    if(in_verbose) cat("YAPSA:::attribute_nucleotide_exchanges::error: ",
                       "in_ALT.field not found. Return NULL.\n")
    return(NULL)
  }
  
  diff <- abs(nchar(in_dat[,REF_ind])-nchar(in_dat[,ALT_ind]))
  
  in_dat$Type <- "dummy"
  for(index in c(1:length(in_dat$CHROM))){
    ref_length <- length(unlist(strsplit(in_dat$REF[index], split="")))
    alt_length <- length(unlist(strsplit(in_dat$ALT[index], split="")))
    if(ref_length > alt_length){
      in_dat$Type[index] <- "Del"
      context_list <- getSequenceContext(position = in_dat$POS[index],
                                         chr = in_dat$CHROM[index],
                                         offsetL= in_offsetL,
                                         offsetR= in_offsetR)
      in_dat$SequenceContext[index] <- context_list$context_string
    }else{
      in_dat$Type[index] <- "Ins"
      context_list <- getSequenceContext(position = in_dat$POS[index],
                                         chr = in_dat$CHROM[index],
                                         offsetL= in_offsetL,
                                         offsetR= in_offsetR)
      in_dat$SequenceContext[index] <- context_list$context_string
    }
  }
  if(length(which(grepl("PID", name_list)))==0){
    in_dat$PID<-"no_PID"
  }
  
  in_dat$Differnce <- diff
  in_dat$Change <- paste0(in_dat[,REF_ind],">",in_dat[,ALT_ind])
  in_dat_return <- in_dat[,c("CHROM", "POS", "REF", "ALT", "PID", "Type",
                             "Differnce", "Change", "SequenceContext")]
  return(in_dat_return)
}

#' Attribution of variant into one onf the 83 INDEL categories
#'
#' Each varaint is categorized into one of the 83 INDEL categories. The
#' classification likewise to Alexandrov et al., 2018
#' (https://www.synapse.org/#!Synapse:syn11726616). The number of 83 features
#' are classefied asfollowed: \enumerate{ \item Deletion of 1 bp C/(G) or T/(A)
#' in a repetitive context. The context is classified into 1, 2, 3, 4, 5 or
#' larger or equal to 6 times the same nucleotide(s). \item Insertion of 1 bp
#' C/(G) or T/(A) in a repetitive context. The context is classified into 0, 1,
#' 2, 3, 4,  or larger or equal to 5 times the same nucleotide(s). \item
#' Deletions of 2bps, 3bps, 4bps or more or equal to 5bps in a repetitive
#' context. Each deletion is classified in a context of 1, 2, 3, 4, 5 or larger
#' or equal to 6 times the same motif. \item Insertion of 2 bps, 3 bps, 4 bps or
#' more or equal to 5 bps in a repetitive context. Each deletion is classified
#' in a context of 0, 1, 2, 3, 4 or larger or equal to 5 times the same motif.
#' \item Microhomology deletion of 2bps, 3bps, 4bps or more or equal to 5 bps in
#' a partly repetitive context. The partly repetitive context is defined by
#' motif length of minus 1 bp, 2 bps, 3 bps, 4 bps or more or equal to 5bps,
#' which is located before and after the break-point junction of the deletion.}
#'
#' @param in_dat_return Data frame constucted form a vcf-like file of a whole
#'   cohort or a single-sample.The first columns are those of a standart vcf
#'   file, followed by an abitrary number of custom or defined columns. One of
#'   these can carry a PID (patient or sample identifyer) and subgroup
#'   information. Furthermore, the columns containing the sequence context and
#'   the absolute length of the INDEL as well as the INDEL type of the variant
#'   can be annotated to the vcf-like df with
#'   \code{\link[YAPSA]{attribute_sequence_contex_indel}}. These columns are
#'   required to enable the constuction of a mutational catalog.
#'
#' @return Data frame with the same dimention as the input data frame  plus an
#'   addional column with the INDEL classification number corrospondig to
#'   Alexandrov et al. 2018.
#'
#' @export
#' @importFrom Biostrings DNAString subseq matchPattern xscat
#'
#' @examples
#' data(GenomeOfNl_raw)
#' GenomeOfNl_context <- attribute_sequence_contex_indel(in_dat =
#' head(GenomeOfNl_raw))
#' GenomeOfNl_classified <- attribution_of_indels(GenomeOfNl_context)
#' GenomeOfNl_classified
#' 
attribution_of_indels <-function(in_dat_return = in_dat_return){ 
  print(paste("INDEL classification of total ", dim(in_dat_return)[1],
              " INDELs This could take a while..."))
  for(indel in c(1:dim(in_dat_return)[1])){
    sequence_string <- DNAString(in_dat_return$SequenceContext[indel])
    if(in_dat_return$Differnce[indel]+11 > 61){
      motive_end_index <- 60
    }else{
      motive_end_index <- in_dat_return$Differnce[indel]+11
    }
    
    if(in_dat_return$Type[indel] == "Del"){
      rest_start_index <- motive_end_index+1
      motive <- sequence_string[12:motive_end_index]
      motive_length <- length(motive)
      rest_string_R <-subseq(sequence_string, start= rest_start_index)
      rest_string_L <-subseq(sequence_string, start=1, end=11)
      match <- matchPattern(motive, rest_string_R)
      
      if(length(width(match))<5){
        rest_string_R_new <- xscat(rest_string_R, DNAString(paste(rep(motive, 5), collapse = '')))
        match <- matchPattern(motive, rest_string_R_new)
      }
      
      if(in_dat_return$Differnce[indel] == 1){
        if(identical(as.character(match),character(0)) | start(match)[1]>1){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 0
            in_dat_return$IndelNumber[indel] <- "DEL_T_1_0"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 0
            in_dat_return$IndelNumber[indel] <- "DEL_C_1_0"  
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length) && 
                 start(match)[4] == 1+(3*motive_length) && 
                 start(match)[5] == 1+(4*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- "5+"
            in_dat_return$IndelNumber[indel] <- "DEL_T_1_5+"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- "5+"
            in_dat_return$IndelNumber[indel] <- "DEL_C_1_5+"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length) && 
                 start(match)[4] == 1+(3*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 4
            in_dat_return$IndelNumber[indel] <- "DEL_T_1_4"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 4
            in_dat_return$IndelNumber[indel] <- "DEL_C_1_4"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 3
            in_dat_return$IndelNumber[indel] <- "DEL_T_1_3"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 3
            in_dat_return$IndelNumber[indel] <- "DEL_C_1_3"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 2
            in_dat_return$IndelNumber[indel] <- "DEL_T_1_2"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 2
            in_dat_return$IndelNumber[indel] <- "DEL_C_1_2"
          }
        }else if(start(match)[1]==1){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 1
            in_dat_return$IndelNumber[indel] <- "DEL_T_1_1"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 1
            in_dat_return$IndelNumber[indel] <- "DEL_C_1_1"
          }
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce[indel] == 2){
        homologR1 <- motive[1:motive_length-1]
        homologL1 <- motive[2:motive_length]
        
        match_homologR1 <- matchPattern(homologR1, rest_string_R)
        match_homologL1 <- matchPattern(homologL1, rest_string_L)
        
        if(((length(as.character(match_homologL1)) != 0 |
             length(as.character(match_homologR1)) != 0) && 
            (start(match)[1]!=1 | identical(as.character(match),
                                            character(0)))) &&
           ((start(match_homologR1)[1]==1 && 
             !is.na(start(match_homologR1)[1]==1))|
            (end(match_homologL1)[1]==11 && 
             !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "1bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_2_1"
        }else if((identical(as.character(match),character(0)) && 
                  is.na(start(match)[1]!=1)) |  start(match)[1]!=1 |
                 identical(as.character(match_homologL1),character(0)) && 
                 identical(as.character(match_homologR1),character(0))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <-"DEL_repeats_2_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) && 
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_2_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) && 
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_2_4"
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_2_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 2
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_2_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_2_1"
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce[indel] == 3){
        homologR1 <- motive[1:motive_length-1]
        homologL1 <- motive[2:motive_length]
        homologR2 <- homologR1[1:length(homologR1)-1]
        homologL2 <- homologL1[2:length(homologL1)]
        
        match_homologR1 <- matchPattern(homologR1, rest_string_R)
        match_homologL1 <- matchPattern(homologL1, rest_string_L)
        match_homologR2 <- matchPattern(homologR2, rest_string_R)
        match_homologL2 <- matchPattern(homologL2, rest_string_L)
        
        
        if(((length(as.character(match_homologL1)) != 0 | 
             length(as.character(match_homologR1)) != 0) && 
            (start(match)[1]!=1 |
             identical(as.character(match),character(0)))) && 
           ((start(match_homologR1)[1]==1 && 
             !is.na(start(match_homologR1)[1]==1))|
            (end(match_homologL1)[1]==11 && 
             !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "2bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_3_2" 
        }else if(((length(as.character(match_homologL2)) != 0 | 
                   length(as.character(match_homologR2)) != 0) && 
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) &&
                 (start(match_homologR2)[1]==1 && 
                  !is.na(start(match_homologR2)[1]==1))|
                 (end(match_homologL2)[1]==11 && 
                  !is.na(end(match_homologL2)[1]==11))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "1bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_3_1" 
        }else if((identical(as.character(match),character(0)) && 
                  is.na(start(match)[1]!=1))| start(match)[1]!=1 | 
                 (identical(as.character(match_homologL1),character(0)) && 
                  identical(as.character(match_homologR1),character(0)) && 
                  identical(as.character(match_homologL2),character(0)) && 
                  identical(as.character(match_homologR2),character(0)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_3_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_3_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_3_4"
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_3_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 2
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_3_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_3_1"
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce[indel] == 4){
        homologR1 <- motive[1:motive_length-1]
        homologL1 <- motive[2:motive_length]
        homologR2 <- homologR1[1:length(homologR1)-1]
        homologL2 <- homologL1[2:length(homologL1)]
        homologR3 <- homologR2[1:length(homologR2)-1]
        homologL3 <- homologL2[2:length(homologL2)]
        
        
        match_homologR1 <- matchPattern(homologR1, rest_string_R)
        match_homologL1 <- matchPattern(homologL1, rest_string_L)
        match_homologR2 <- matchPattern(homologR2, rest_string_R)
        match_homologL2 <- matchPattern(homologL2, rest_string_L)
        match_homologR3 <- matchPattern(homologR3, rest_string_R)
        match_homologL3 <- matchPattern(homologL3, rest_string_L)
        
        
        if(((length(as.character(match_homologL1)) != 0 |
             length(as.character(match_homologR1)) != 0) &&
            (start(match)[1]!=1 | 
             identical(as.character(match),character(0)))) && 
           ((start(match_homologR1)[1]==1 && 
             !is.na(start(match_homologR1)[1]==1))|
            (end(match_homologL1)[1]==11 && 
             !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "3bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_4_3"
        }else if(((length(as.character(match_homologL2)) != 0 | 
                   length(as.character(match_homologR2)) != 0) &&
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) && 
                 ((start(match_homologR2)[1]==1 &&
                   !is.na(start(match_homologR2)[1]==1))|
                  (end(match_homologL2)[1]==11 && 
                   !is.na(end(match_homologL2)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "2bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_4_2"
        }else if(((length(as.character(match_homologL3)) != 0 | 
                   length(as.character(match_homologR3)) != 0) &&
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) && 
                 ((start(match_homologR3)[1]==1 && 
                   !is.na(start(match_homologR3)[1]==1))|
                  (end(match_homologL3)[1]==11 && 
                   !is.na(end(match_homologL3)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "1bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_4_1"
        }else if((identical(as.character(match),character(0)) && 
                  is.na(start(match)[1]!=1)) |  start(match)[1]!=1 |
                 (identical(as.character(match_homologL1),character(0)) && 
                  identical(as.character(match_homologR1),character(0)) && 
                  identical(as.character(match_homologL2),character(0)) &&
                  identical(as.character(match_homologR2),character(0)) &&
                  identical(as.character(match_homologL3),character(0)) && 
                  identical(as.character(match_homologR3),character(0)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_4_0" 
        }else if(length(match) >= 5 && 
                 (start(match)[1] == 1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) && 
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_4_5+" 
        }else if(length(match) >= 4 &&
                 (start(match)[1] == 1 &&
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_4_4" 
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_4_3" 
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 2
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_4_2" 
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_4_1" 
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce[indel] >= 5){
        
        motive_first_5bp <- motive[1:5]
        start_last_5bp <- (motive_length-5)+1
        motive_last_5bp <- motive[start_last_5bp:motive_length]
        
        homologR1 <- motive_first_5bp
        homologL1 <- motive_last_5bp
        homologR2 <- homologR1[1:length(homologR1)-1]
        homologL2 <- homologL1[2:length(homologL1)]
        homologR3 <- homologR2[1:length(homologR2)-1]
        homologL3 <- homologL2[2:length(homologL2)]
        homologR4 <- homologR3[1:length(homologR3)-1]
        homologL4 <- homologL3[2:length(homologL3)]
        homologR5 <- homologR4[1:length(homologR4)-1]
        homologL5 <- homologL4[2:length(homologL4)]
        
        
        match_homologR1 <- matchPattern(homologR1, rest_string_R)
        match_homologL1 <- matchPattern(homologL1, rest_string_L)
        match_homologR2 <- matchPattern(homologR2, rest_string_R)
        match_homologL2 <- matchPattern(homologL2, rest_string_L)
        match_homologR3 <- matchPattern(homologR3, rest_string_R)
        match_homologL3 <- matchPattern(homologL3, rest_string_L)
        match_homologR4 <- matchPattern(homologR4, rest_string_R)
        match_homologL4 <- matchPattern(homologL4, rest_string_L)
        match_homologR5 <- matchPattern(homologR5, rest_string_R)
        match_homologL5 <- matchPattern(homologL5, rest_string_L)
        
        if(length(match) >= 5 && 
           (start(match)[1]==1 &&
            start(match)[2] == 1+motive_length &&
            start(match)[3] == 1+(2*motive_length) &&
            start(match)[4] == 1+(3*motive_length) &&
            start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_5+_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) && 
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_5+_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_5+_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_5+_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_5+_1"
        }else if(((length(as.character(match_homologL1)) != 0 | 
                   length(as.character(match_homologR1)) != 0) && 
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) &&
                 ((start(match_homologR1)[1]==1 && 
                   !is.na(start(match_homologR1)[1]==1))|
                  (end(match_homologL1)[1]==11 && 
                   !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_5+_5+"
        }else if(((length(as.character(match_homologL2)) != 0 | 
                   length(as.character(match_homologR2)) != 0) && 
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) &&
                 ((start(match_homologR2)[1]==1 && 
                   !is.na(start(match_homologR2)[1]==1))|
                  (end(match_homologL2)[1]==11 && 
                   !is.na(end(match_homologL2)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "4bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_5+_4"
        }else if(((length(as.character(match_homologL3)) != 0 |
                   length(as.character(match_homologR3)) != 0) &&
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0))))&&
                 ((start(match_homologR3)[1]==1 &&
                   !is.na(start(match_homologR3)[1]==1))|
                  (end(match_homologL3)[1]==11 &&
                   !is.na(end(match_homologL3)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "3bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_5+_3"
        }else if(((length(as.character(match_homologL4)) != 0 |
                   length(as.character(match_homologR4)) != 0) && 
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0))))&&
                 ((start(match_homologR4)[1]==1 && 
                   !is.na(start(match_homologR4)[1]==1))|
                  (end(match_homologL4)[1]==11 &&
                   !is.na(end(match_homologL4)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "2bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_5+_2"
        }else if(((length(as.character(match_homologL5)) != 0 |
                   length(as.character(match_homologR5)) != 0) &&
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) &&
                 ((start(match_homologR5)[1]==1 && 
                   !is.na(start(match_homologR5)[1]==1))|
                  (end(match_homologL5)[1]==11 &&
                   !is.na(end(match_homologL5)[1]==11)))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "1bp"
          in_dat_return$IndelNumber[indel] <- "DEL_MH_5+_1"
        }else if(((identical(as.character(match),character(0)) && 
                   is.na(start(match)[1]!=1)) |  start(match)[1]!=1 | 
                  identical(as.character(match_homologL1),character(0)) &&
                  identical(as.character(match_homologR1),character(0)) &&
                  identical(as.character(match_homologL2),character(0)) && 
                  identical(as.character(match_homologR2),character(0)) && 
                  identical(as.character(match_homologL3),character(0)) && 
                  identical(as.character(match_homologR3),character(0)) &&
                  identical(as.character(match_homologL4),character(0)) && 
                  identical(as.character(match_homologR4),character(0)) &&
                 identical(as.character(match_homologL5),character(0)) && 
                 identical(as.character(match_homologR5),character(0)))
        ){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <- "DEL_repeats_5+_0"
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }
    }else{
      rest_start_index <- 12
      motive <- substr(in_dat_return$ALT[indel], start = 2, stop = nchar(in_dat_return$ALT[indel]))
      motive_length <- length(motive)
      rest_string_R <-subseq(sequence_string, start = rest_start_index)
      rest_string_L <-subseq(sequence_string, start=1, end=11)
      match <- matchPattern(motive, rest_string_R)
      
      if(length(width(match))<5){
        rest_string_R_new <- xscat(rest_string_R, DNAString(paste(rep(motive, 5), collapse = '')))
        match <- matchPattern(motive, rest_string_R_new)
      }
      
      if(in_dat_return$Differnce[indel] == 1){
        if(start(match)[1]==1 && 
           start(match)[2] == 1+motive_length && 
           start(match)[3] == 1+(2*motive_length) && 
           start(match)[4] == 1+(3*motive_length) && 
           start(match)[5] == 1+(4*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- "5+"
            in_dat_return$IndelNumber[indel] <- "INS_T_1_5+"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- "5+"
            in_dat_return$IndelNumber[indel] <- "INS_C_1_5+"
          }
        }else if(start(match)[1]==1 && 
                 start(match)[2] == 1+motive_length &&
                 start(match)[3] == 1+(2*motive_length) && 
                 start(match)[4] == 1+(3*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 4
            in_dat_return$IndelNumber[indel] <- "INS_T_1_4"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 4
            in_dat_return$IndelNumber[indel] <- "INS_C_1_4"
          }
        }else if(start(match)[1]==1 && 
                 start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 3
            in_dat_return$IndelNumber[indel] <- "INS_T_1_3"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 3
            in_dat_return$IndelNumber[indel] <- "INS_C_1_3"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 2
            in_dat_return$IndelNumber[indel] <- "INS_T_1_2"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 2
            in_dat_return$IndelNumber[indel] <- "INS_C_1_2"
          }
        }else if(start(match)[1]==1){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 1
            in_dat_return$IndelNumber[indel] <- "INS_T_1_1"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 1
            in_dat_return$IndelNumber[indel] <- "INS_C_1_1"
          }
        }else{
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv[indel] <- "T:A"
            in_dat_return$RepeatSize[indel] <- 0
            in_dat_return$IndelNumber[indel] <- "INS_T_1_0"
          }else{
            in_dat_return$Motiv[indel] <- "C:G"
            in_dat_return$RepeatSize[indel] <- 0
            in_dat_return$IndelNumber[indel] <- "INS_C_1_0"  
          }
        }
      }else if(in_dat_return$Differnce[indel] == 2){
        if (identical(as.character(match),character(0)) | start(match)[1]!=1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <- "INS_repeats_2_0"
        }else if(length(match) >= 5 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length) &&
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "INS_repeats_2_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "INS_repeats_2_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "INS_repeats_2_3"
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 2
          in_dat_return$IndelNumber[indel] <- "INS_repeats_2_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "INS_repeats_2_1"
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce[indel] == 3){
        if (identical(as.character(match),character(0))| start(match)[1]!=1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <- "INS_repeats_3_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "INS_repeats_3_5+"
        }else if (length(match) >= 4 && 
                  (start(match)[1]==1 && 
                   start(match)[2] == 1+motive_length && 
                   start(match)[3] == 1+(2*motive_length) &&
                   start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "INS_repeats_3_4"
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "INS_repeats_3_3"
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 2
          in_dat_return$IndelNumber[indel] <- "INS_repeats_3_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "INS_repeats_3_1"
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce[indel] == 4){
        if (identical(as.character(match),character(0)) | start(match)[1]!=1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <- "INS_repeats_4_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length) &&
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "INS_repeats_4_5+"
        }else if(length(match) >= 4 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "INS_repeats_4_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "INS_repeats_4_3"
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 2
          in_dat_return$IndelNumber[indel] <- "INS_repeats_4_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "INS_repeats_4_1"
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce[indel] >= 5){
        if (identical(as.character(match),character(0))| start(match)[1]!=1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 0
          in_dat_return$IndelNumber[indel] <- "INS_repeats_5+_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- "5+"
          in_dat_return$IndelNumber[indel] <- "INS_repeats_5+_5+"
        }else if(length(match) >= 4 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 4
          in_dat_return$IndelNumber[indel] <- "INS_repeats_5+_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 3
          in_dat_return$IndelNumber[indel] <- "INS_repeats_5+_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 2
          in_dat_return$IndelNumber[indel] <- "INS_repeats_5+_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv[indel] <- as.character(motive)
          in_dat_return$RepeatSize[indel] <- 1
          in_dat_return$IndelNumber[indel] <- "INS_repeats_5+_1"
        }else{
          in_dat_return$Motiv[indel] <- "NA"
          in_dat_return$RepeatSize[indel] <- "NA"
          in_dat_return$IndelNumber[indel] <- "NOT DEFINED"
        }
      }
    }  
  }   
  return(in_dat_return)
}  

#' Plot the spectra of nucleotide exchanges of INDELs
#'
#' Plots the spectra of nucelotides in their triplet contexts. If several
#' columns are present in the input data frame, the spectra are ploted for every
#' column seperatly. The function is only suitable for a INDEL spectra and for
#' SNV representation the funtion \code{\link[YAPSA]{plotExchangeSpectra}}
#' should be used.
#'
#' @param in_catalogue_df Numerical data frame encoding the exchange spectra to
#'   be displayed, either a mutational catalogue \code{V} or a signatures matrix
#'   \code{W}
#'
#' @param in_colour_vector Specifies the colours of the INDELs if non-null
#' @param in_show_indel Whether or not to show the INDEL names on the x-axis
#' @param in_show_axis_title Whether or not to show the name of the y-axis
#' @param in_scales Argument passed on to \code{\link[ggplot2]{facet_grid}}
#' @param in_refLine If non-null, value on the y-axis at which a horizontal line
#'   is to be drawn
#' @param in_refAlpha Transparency of the horizontal line if it is to be drawn
#' @param in_background Option to provide a background theme, e.g.
#'   \code{\link[ggplot2]{theme_grey}}

#'
#' @return The generated barplot - a ggplot2 plot
#'
#' @export
#' @importFrom reshape2 melt
#' @import ggplot2
#' @examples
#' data(sigs_pcawg)
#' plotExchangeSpectra_indel(PCAWG_SP_ID_sigs_df[,c(6,8)])
#' 
plotExchangeSpectra_indel <- function (in_catalogue_df, 
                                       in_colour_vector = NULL,
                                       in_show_indel = FALSE, 
                                       in_show_axis_title = FALSE, 
                                       in_scales = "free_x",
                                       in_refLine = NULL, 
                                       in_refAlpha = 0.5, 
                                       in_background = NULL){
  .e <- environment()
  in_catalogue_df$indel_exchange <- rownames(in_catalogue_df)
  in_catalogue_df$nuc_exchange <- c(rep("DEL_C", 6),rep("DEL_T",6),
                                    rep("INS_C",6),rep("INS_T",6),
                                    rep("DEL_2", 6),rep("DEL_3", 6),
                                    rep("DEL_4", 6),rep("DEL_5+", 6),
                                    rep("INS_2", 6),rep("INS_3", 6),
                                    rep("INS_4", 6),rep("INS_5+", 6),
                                    rep("MH_2", 1),rep("MH_3", 2),
                                    rep("MH_4", 3),rep("MH_5+", 5))
  in_catalogue_df$repetion <- c(rep(c(1,2,3,4,5,"6+"), 2), 
                                rep(c(0,1,2,3,4,"5+"), 2), 
                                rep(c(1,2,3,4,5,"6+"), 4),
                                rep(c(0,1,2,3,4,"5+"), 4), 
                                1,1,2,1,2,3,1,2,3,4,"5+")
  catalogue_df_melt <- melt(in_catalogue_df, id.vars = c("nuc_exchange", 
                                                         "repetion", 
                                                         "indel_exchange"))
  catalogue_df_melt$nuc_exchange <- factor(catalogue_df_melt$nuc_exchange, 
                                           levels = unique(catalogue_df_melt$nuc_exchange))
  
  my_palette <- c("tan1", "darkorange1", "darkseagreen3", "forestgreen", 
                  "rosybrown1", "salmon", "red3","darkred",
                  "#C9DFF3", "skyblue2", "steelblue3", "steelblue4",
                  "#D4D9FF", "#BFC1E6","#8E8EBE", "mediumpurple4")
  names(my_palette) <- c("DEL_C", "DEL_T", "INS_C", "INS_T",
                         paste0("DEL_", c(2:4, "5+")),
                         paste0("INS_", c(2:4, "5+")),
                         paste0("MH_", c(2:4, "5+")))
  if (!is.null(in_colour_vector)) 
    my_palette <- in_colour_vector
  
  p <- ggplot(catalogue_df_melt, 
              environment = .e) + geom_bar(aes_string(x = "repetion", 
                                                      y = "value",
                                                      fill = "nuc_exchange"),
                                           stat = "identity") + 
    scale_fill_manual(name = "exchange", values = my_palette) + 
    facet_grid(variable ~ nuc_exchange, scales = in_scales, space = in_scales)+
    theme(strip.text.x = element_text(size=7))
  #theme(strip.background = element_rect(fill = my_palette))
  if (!is.null(in_background)) {
    p <- p + in_background
  }
  if (is.numeric(in_refLine)) {
    p <- p + geom_hline(yintercept = in_refLine, alpha = in_refAlpha)
  }
  p1 <- p + theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank())
  if (in_show_indel) {
    p1 <- p + theme(axis.title.x = element_blank(), 
                    axis.text.x = element_text(angle = 0, 
                                               vjust = 0.5), axis.ticks.x = element_blank())
  }
  if (!in_show_axis_title) 
    p1 <- p1 + theme(axis.title.y = element_blank())
  return(p1)
}  

#' Create a Mutational catalog from a data frame
#'
#' This function creates a mutational catalog from a data frame. It requires the
#' returend data frame optainted with
#' \code{\link[YAPSA]{attribution_of_indels}}.
#'
#' @param in_df A data frame constucted from a vcf-like file of a whole cohort
#'   or single-sample. The first coloums are those of a standart vcf file,
#'   followed by an arbitrary number of customs or used defined columns. One if
#'   these can carry a PID (patient or sample identefyier) and the subgroup
#'   information. Additionaly to consuct the the mutational catalog each variant
#'   needs to be characterize into one of the 83 INDEL feature classes, which
#'   can be perfomed with \code{\link[YAPSA]{attribution_of_indels}}
#' @param in_signatures_df A numeric data frame \code{W} with \code{n} rows and
#'   \code{l} columns, \code{n} being the number of features and \code{l} being
#'   the number of signatures.Data frame containing INDEL signatures which
#'   should be used to create the mutational cataolog \code{V}.
#'
#' @return A count dataframe, the mutational catalog \code{V} with rownames
#'   indicating the INDELs and colnames having the PIDs
#'
#' @export
#' @importFrom reshape2 melt
#' @importFrom Biostrings DNAString subseq matchPattern xscat
#'
#' @examples
#' data(GenomeOfNl_raw)
#' data(sigs_pcawg)
#' GenomeOfNl_context <- attribute_sequence_contex_indel(in_dat =
#'  head(GenomeOfNl_raw))
#' GenomeOfNl_classified <- attribution_of_indels(GenomeOfNl_context)
#' GenomeOfNl_mut_cat <- create_indel_mut_cat_from_df(GenomeOfNl_classified,
#'  in_signatures_df=PCAWG_SP_ID_sigs_df)
#' 
create_indel_mut_cat_from_df <- function(in_df, in_signatures_df = NULL){
  PID <- unique(in_df$PID)
  
  if(is.null(in_signatures_df)){
    channels <- c('DEL_C_1_0', 'DEL_C_1_1', 'DEL_C_1_2', 'DEL_C_1_3', 'DEL_C_1_4', 'DEL_C_1_5+', 'DEL_T_1_0', 'DEL_T_1_1', 'DEL_T_1_2', 'DEL_T_1_3', 'DEL_T_1_4', 'DEL_T_1_5+',
                  'INS_C_1_0', 'INS_C_1_1', 'INS_C_1_2', 'INS_C_1_3', 'INS_C_1_4', 'INS_C_1_5+', 'INS_T_1_0', 'INS_T_1_1', 'INS_T_1_2', 'INS_T_1_3', 'INS_T_1_4', 'INS_T_1_5+',
                  'DEL_repeats_2_0', 'DEL_repeats_2_1', 'DEL_repeats_2_2', 'DEL_repeats_2_3', 'DEL_repeats_2_4', 'DEL_repeats_2_5+', 'DEL_repeats_3_0', 'DEL_repeats_3_1', 
                  'DEL_repeats_3_2', 'DEL_repeats_3_3', 'DEL_repeats_3_4', 'DEL_repeats_3_5+', 'DEL_repeats_4_0', 'DEL_repeats_4_1', 'DEL_repeats_4_2', 'DEL_repeats_4_3',
                  'DEL_repeats_4_4', 'DEL_repeats_4_5+', 'DEL_repeats_5+_0', 'DEL_repeats_5+_1', 'DEL_repeats_5+_2', 'DEL_repeats_5+_3', 'DEL_repeats_5+_4', 'DEL_repeats_5+_5+',
                  'INS_repeats_2_0', 'INS_repeats_2_1', 'INS_repeats_2_2', 'INS_repeats_2_3', 'INS_repeats_2_4', 'INS_repeats_2_5+', 'INS_repeats_3_0', 'INS_repeats_3_1',
                  'INS_repeats_3_2', 'INS_repeats_3_3', 'INS_repeats_3_4', 'INS_repeats_3_5+', 'INS_repeats_4_0', 'INS_repeats_4_1', 'INS_repeats_4_2', 'INS_repeats_4_3',
                  'INS_repeats_4_4', 'INS_repeats_4_5+', 'INS_repeats_5+_0', 'INS_repeats_5+_1', 'INS_repeats_5+_2', 'INS_repeats_5+_3', 'INS_repeats_5+_4', 'INS_repeats_5+_5+',
                  'DEL_MH_2_1', 'DEL_MH_3_1', 'DEL_MH_3_2', 'DEL_MH_4_1', 'DEL_MH_4_2', 'DEL_MH_4_3', 'DEL_MH_5+_1', 'DEL_MH_5+_2', 'DEL_MH_5+_3', 'DEL_MH_5+_4', 'DEL_MH_5+_5+')
    count_df <- data.frame(matrix(data = 0, nrow = length(channels), 
                                  ncol=0), row.names = channels)
  } else {
    count_df <- data.frame(matrix(data = 0, nrow = dim(in_signatures_df)[1], 
                                  ncol=0), row.names = rownames(in_signatures_df))
  }
  
  for(pid in PID){
    in_subset_df <- subset(in_df, in_df$PID %in% list(pid), drop =TRUE)
    count_per_pid <- data.frame(table(in_subset_df$IndelNumber))
    rownames(count_per_pid) <- count_per_pid$Var1
    count_per_pid$Var1 <- NULL
    colnames(count_per_pid) <- pid
    count_df <- setNames(
      cbind(count_df, 
            count_per_pid[, pid][match(rownames(count_df),
                                       rownames(count_per_pid))]),
      c(colnames(count_df),pid))
  }
  count_df[is.na(count_df)] <- 0
  return(count_df)
}

#' Wrapper to create an INDEL mutational catalog from a vlf-like data frame
#'
#' From data frame constucted from a vcf-file file the function
#' \code{\link[YAPSA]{create_indel_mutation_catalogue_from_df}} creates a
#' mutational catalog V by squencially applying the
#' \code{\link[YAPSA]{attribute_sequence_contex_indel}},
#' \code{\link[YAPSA]{attribute_sequence_contex_indel}} and then
#' \code{\link[YAPSA]{attribution_of_indels}}. The runtime of the function is
#' about 1 sec per 6 variants as sequence context as well as INDEL
#' calssification are timeconsuming to compute (optimization ongoing)
#'
#' @param in_dat A data frame constructed from a vcf-like file of a whole cohort
#'   or single-sample. The first columns are those of a standard vcf file
#'   (\code{CHROM}, \code{POS}, \code{REF} and \code{ALT}), followed by an
#'   arbitrary number of custom or used defined columns. One of these can carry
#'   a PID (patient or sample identifyier) and one can carry subgroup
#'   information.
#' @param in_signature_df A numeric data frame \code{W} with \code{n} rows and
#'   \code{l} columns, n being the number of features and l being the number od
#'   signatures. Data frame containing INDEL signatures which should be used to
#'   create the mutational cataologe \code{V}.
#' @param in_REF.field String indicating which column of \code{in_dat} carries
#'   the reference base if dealing with data frames
#' @param in_ALT.field String indicating which column of \code{in_dat} carries
#'   the variant base if dealing with data frames
#' @param in_verbose Verbose if \code{in_verbose=1}
#'
#'
#' @return A dataframe in the format of a mutational catalog \code{V}, which can
#' be used for \code{\link[YAPSA]{LCD}} analysis
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom Biostrings getSeq DNAString subseq matchPattern xscat
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export
#'
#' @examples
#' data(sigs_pcawg)
#' data(GenomeOfNl_raw)
#' temp_df <- translate_to_hg19(GenomeOfNl_raw[1:200,],"CHROM")
#' temp_df$PID <- sample(c("PID1","PID2","PID3","PID4","PID5"),200,replace=TRUE)
#' temp <- create_indel_mutation_catalogue_from_df(in_dat = temp_df,
#'    in_signature_df = PCAWG_SP_ID_sigs_df,
#'    in_REF.field = "REF",
#'    in_ALT.field = "ALT",
#'    in_verbose = FALSE)
#' dim(temp)
#' head(temp)
#' 
create_indel_mutation_catalogue_from_df <- function(in_dat, 
                                                    in_signature_df = NULL, 
                                                    in_REF.field="REF",
                                                    in_ALT.field="ALT",
                                                    in_verbose = FALSE){
  indel_context_atrribution <- attribute_sequence_contex_indel(
    in_dat, 
    in_REF.field=in_REF.field,
    in_ALT.field=in_ALT.field, 
    in_verbose = FALSE)
  indel_classification_attribution <- attribution_of_indels(
    in_dat_return = indel_context_atrribution)
  indel_mut_cat_creation <- create_indel_mut_cat_from_df(
    indel_classification_attribution, in_signature_df)
  return(indel_mut_cat_creation)
}


#' Wrapper to compute confidence intervals for SNV and INDEL signatures of a
#' cohort or single-sample
#'
#' Wrapper function around \code{\link[YAPSA]{confIntExp}}, which is applies to
#' every signature or sample pair in a cohort. The extracted lower bound of the
#' confidence intervals are added to the input data which is reodered and melted
#' in order to prepare for visualization with ggplot2. The calculates of
#' confidence intervals is based on a profiling likelihood algorithm and the
#' wrapper calculates the data for the exposure contubution identefied with SNV
#' and INDEL signature decompositions and application of the following cutoffs:
#' \enumerate{ \item \code{CosmicValid_absCutoffVector} \item
#' \code{CosmicValid_normCutoffVector} \item \code{CosmicArtif_absCutoffVector}
#' \item \code{CosmicArtif_normCutoffVector} \item
#' \code{PCAWGValidSNV_absCutoffVector} \item
#' \code{PCAWGValidID_absCutoffVector} } The function makes use of differnet
#' YAPSA functions. For each of the above stated cutoff vectors a per PID
#' decompostion of the SNV and INDEL catalog is calulated respectivly using
#' \code{\link[YAPSA]{LCD_complex_cutoff_perPID}}. In a next step,
#' \code{\link[YAPSA]{variateExp}} wich is a wrapper around
#' \code{\link[YAPSA]{confIntExp}} to compute confidence intervals for a cohort
#' is used. A dataframe is returend with the upper and lower bounds of the
#' confidence intervals. In a last step
#' \code{\link[YAPSA]{plotExposuresConfidence_indel}} to plot the exposures to
#' extracted signatures including confidence intervals computed with e.g. by
#' \code{\link[YAPSA]{variateExp}}.
#'
#'
#' @param in_current_indel_df A INDEL mutational catalog. Mutational catalog can
#'   be constucted with
#'   \code{\link[YAPSA]{create_indel_mutation_catalogue_from_df}}
#' @param in_current_snv_df A SNV mutational catalog. Mutational catalog can be
#'   constuced with \code{\link[YAPSA]{create_mutation_catalogue_from_df}}
#'
#' @return A list is returned containing 12 objects. For each cutoff data frame
#'   two corrosponding object are present. First, the \code{p} gtable object
#'   which can be used for gaphically visualization, and second a dataframe
#'   containing  the corrosponding upper and lower bounds of the confidence
#'   intervals.
#'
#' @export
#'
#' @examples
#' data("GenomeOfNl_MutCat")
#' 
confidence_indel_calulation <- function(in_current_indel_df, 
                                        in_current_snv_df){
  data(sigs_pcawg)
  data(cutoffs_pcawg)
  data(sigs)
  data(cutoffs)
  
  CosmicValid_absCutoffVector <- cutoffCosmicValid_abs_df[6, ]
  CosmicValid_normCutoffVector <- cutoffCosmicValid_rel_df[6, ]
  CosmicArtif_absCutoffVector <- cutoffCosmicArtif_abs_df[6, ]
  CosmicArtif_normCutoffVector <- cutoffCosmicArtif_rel_df[6, ]
  
  PCAWGValidSNV_absCutoffVector <- cutoffPCAWG_SBS_WGSWES_realPid_df[13, ]
  PCAWGValidID_absCutoffVector <- cutoffPCAWG_ID_WGS_Pid_df[3, ]
  
  current_id_catalogue_df <- data.frame(in_current_indel_df)
  current_snv_catalogue_df <- data.frame(in_current_snv_df)
  
  
  PCAWGValidSNV_abs_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df = current_snv_catalogue_df,
    in_signatures_df = PCAWG_SP_SBS_sigs_Real_df,
    in_cutoff_vector = PCAWGValidSNV_absCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind_df = PCAWG_SP_SBS_sigInd_Real_df)
  
  PCAWGValidID_abs_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df = current_id_catalogue_df,
    in_signatures_df = PCAWG_SP_ID_sigs_df,
    in_cutoff_vector = PCAWGValidID_absCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind_df = PCAWG_SP_ID_sigInd_df)
  
  CosmicValid_abs_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df = current_snv_catalogue_df,
    in_signatures_df = AlexCosmicValid_sig_df,
    in_cutoff_vector = CosmicValid_absCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind_df = AlexCosmicValid_sigInd_df)
  
  CosmicValid_norm_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df = current_snv_catalogue_df,
    in_signatures_df = AlexCosmicValid_sig_df,
    in_cutoff_vector = CosmicValid_normCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind_df = AlexCosmicValid_sigInd_df)
  
  CosmicArtif_abs_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df = current_snv_catalogue_df,
    in_signatures_df = AlexCosmicArtif_sig_df,
    in_cutoff_vector = CosmicArtif_absCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind_df = AlexCosmicArtif_sigInd_df)
  
  CosmicArtif_norm_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df = current_snv_catalogue_df,
    in_signatures_df = AlexCosmicArtif_sig_df,
    in_cutoff_vector = CosmicArtif_normCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind_df = AlexCosmicArtif_sigInd_df)
  
  
  
  
  exposures_list <- list(PCAWGValidSNV_abs =PCAWGValidSNV_abs_LCDlist$exposures,
                         PCAWGValidID_abs = PCAWGValidID_abs_LCDlist$exposures,
                         CosmicValid_abs = CosmicValid_abs_LCDlist$exposures,
                         CosmicValid_norm = CosmicValid_norm_LCDlist$exposures,
                         CosmicArtif_abs = CosmicArtif_abs_LCDlist$exposures,
                         CosmicArtif_norm = CosmicArtif_norm_LCDlist$exposures)
  
  name_list <- names(exposures_list)
  exposures_list <- lapply(seq_along(exposures_list),FUN=function(current_ind){
    current_df <- exposures_list[[current_ind]]
    names(current_df) <- name_list[current_ind]
    return(current_df)
  })
  
  
  names(exposures_list) <- name_list
  
  exposures_list <- lapply(exposures_list, function(list){
    colnames(list) <- colnames(current_id_catalogue_df)
    return(list)
  })
  
  subgroups_df <- data.frame(
    PID = colnames(exposures_list[[1]]), subgroup = "test", sum = 1, 
    compl_sum = 1, index = seq_along(exposures_list[[1]]), col = "#FF0000FF")
  
  
  
  complete_df_list <- lapply(name_list, function(current_condition){
    print(current_condition)
    if(grepl("CosmicValid_abs|CosmicValid_norm|CosmicArtif_abs|CosmicArtif_norm",
             current_condition)){
      current_exposures_df <- exposures_list[[current_condition]]
      current_sigs <- rownames(current_exposures_df)
      current_signatures_df <- AlexCosmicArtif_sig_df[, current_sigs,
                                                      drop=FALSE]
      suppressWarnings(variateExp(
        in_catalogue_df = current_snv_catalogue_df,
        in_sig_df = current_signatures_df,
        in_exposures_df = current_exposures_df,
        in_sigLevel = 0.025, in_delta = 0.4))
    }else if(grepl("PCAWGValidSNV_abs|PCAWGValidSNV_norm", current_condition)){
      current_exposures_df <- exposures_list[[current_condition]]
      current_sigs <- rownames(current_exposures_df)
      current_signatures_df <- PCAWG_SP_SBS_sigs_Real_df[, current_sigs, 
                                                         drop=FALSE]
      suppressWarnings(variateExp(
        in_catalogue_df = current_snv_catalogue_df,
        in_sig_df = current_signatures_df,
        in_exposures_df = current_exposures_df,
        in_sigLevel = 0.025, in_delta = 0.4))
    }else if(grepl("PCAWGValidID_abs", current_condition)){
      current_exposures_df <- exposures_list[[current_condition]]
      current_sigs <- rownames(current_exposures_df)
      current_signatures_df <- PCAWG_SP_ID_sigs_df[, current_sigs, 
                                                   drop=FALSE]
      suppressWarnings(variateExp(
        in_catalogue_df = current_id_catalogue_df,
        in_sig_df = current_signatures_df,
        in_exposures_df = current_exposures_df,
        in_sigLevel = 0.025, in_delta = 0.4))
    }else{
      print("WARNING: There is no match between names of exposure list and the
            current condition")
    }
    })
  
  complete_df_list <- lapply(complete_df_list, function(current_complete_df){
    current_complete_df$norm_exposure <- current_complete_df$exposure / 
      current_complete_df[which(current_complete_df$sig == "total"), "exposure"]
    current_complete_df$norm_lower <- current_complete_df$norm_exposure * 
      current_complete_df$relLower
    current_complete_df$norm_upper <- current_complete_df$norm_exposure * 
      current_complete_df$relUpper
    return(current_complete_df)
  })
  
  names(complete_df_list)<-names(exposures_list)
  
  one <- plotExposuresConfidence_indel(complete_df_list$PCAWGValidSNV_abs, 
                                       subgroups_df,
                                       PCAWG_SP_SBS_sigInd_Real_df)
  two <- plotExposuresConfidence_indel(complete_df_list$PCAWGValidID_abs, 
                                       subgroups_df,
                                       PCAWG_SP_ID_sigInd_df)
  three <- plotExposuresConfidence_indel(complete_df_list$CosmicValid_abs, 
                                         subgroups_df,
                                         AlexCosmicArtif_sigInd_df)
  four <- plotExposuresConfidence_indel(complete_df_list$CosmicArtif_abs, 
                                        subgroups_df,
                                        AlexCosmicArtif_sigInd_df)
  five <- plotExposuresConfidence_indel(complete_df_list$CosmicValid_norm, 
                                        subgroups_df,
                                        AlexCosmicArtif_sigInd_df)
  six <- plotExposuresConfidence_indel(complete_df_list$CosmicArtif_norm,
                                       subgroups_df,
                                       AlexCosmicArtif_sigInd_df)
  
  
  return(list(p_complete_PCAWG_SNV = one,
              p_complete_PCAWG_ID = two,
              p_complete_COSMIC_valid_abs = three, 
              p_complete_COSMIC_artif_abs = four,
              p_complete_COSMIC_valid_norm = five, 
              p_complete_COSMIC_artif_norm = six,
              complete_PCAWG_SNV=complete_df_list$PCAWGValidSNV_abs,
              complete_PCAWG_ID=complete_df_list$PCAWGValidID_abs,
              complete_COSMIC_valid_abs = complete_df_list$CosmicValid_abs,
              complete_COSMIC_artif_abs = complete_df_list$CosmicArtif_abs,
              complete_COSMIC_valid_norm = complete_df_list$CosmicValid_norm,
              complete_COSMIC_artif_norm = complete_df_list$CosmicArtif_norm))
  }



#' Wrapper to compute confidence intervals for only INDEL signatures.
#'
#' Wrapper function around \code{\link[YAPSA]{confIntExp}}, which is applies to
#' every signature or sample pair in a cohort. The extracted lower bound of the
#' confidence intervals are added to the input data which is reodered and melted
#' in order to prepare for visualization with ggplot2. The calculates of
#' confidence intervals is based on a profiling likelihood algorithm and the
#' wrapper calculates the data for the exposure contubution identefied with
#' INDEL singature decomposition and the usage of
#' \code{PCAWGValidID_absCutoffVector} data frame.
#'
#' The function makes use of differnet YAPSA functions. For each of the above
#' stated cutoff vectors a per PID decompostion of the SNV and INDEL catalog is
#' calulated respectivly using \code{\link[YAPSA]{LCD_complex_cutoff_perPID}}.
#' In a next step, \code{\link[YAPSA]{variateExp}} which is a wrapper around
#' \code{\link[YAPSA]{confIntExp}} to compute confidenceintervals for a cohort
#' is used. A dataframe is returend with the upper and lower bounds of the
#' confidence intervals. In a last step
#' \code{\link[YAPSA]{plotExposuresConfidence_indel}} to plot the exposures to
#' extracted signatures including confidence intervals computed with e.g. by
#' \code{\link[YAPSA]{variateExp}}.
#'
#'
#' @param in_current_indel_df A INDEL mutational catalog. Mutational catalog can
#'   be constucted with \code{\link{create_indel_mutation_catalogue_from_df}}
#'
#' @return A list is returned containing two object. First, the \code{p} gtable
#' object which can be used for gaphically visualization, and second a dataframe
#' containing the corrosponding upper and lower bounds of the confidence
#' intervals.
#'
#' @export
#'
#' @examples
#' data("GenomeOfNl_MutCat")
#' temp_list <- confidence_indel_only_calulation(
#'                         in_current_indel_df=MutCat_indel_df)
#' plot(temp_list$p_complete_PCAWG_ID)
#' head(temp_list$complete_PCAWG_ID)
#' 
confidence_indel_only_calulation <- function(in_current_indel_df){
  
  data(sigs_pcawg)
  data(cutoffs_pcawg)
  #devtools:::use_data(sigs_pcawg, overwrite = FALSE)
  #devtools:::use_data(cutoffs_pcawg, overwrite = FALSE)
  
  PCAWGValidID_absCutoffVector <- cutoffPCAWG_ID_WGS_Pid_df[3, ]
  
  current_id_catalogue_df <- data.frame(in_current_indel_df, check.names=FALSE)
  
  PCAWGValidID_abs_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df = current_id_catalogue_df,
    in_signatures_df = PCAWG_SP_ID_sigs_df,
    in_cutoff_vector = PCAWGValidID_absCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind_df = PCAWG_SP_ID_sigInd_df)
  
  
  exposures_list <- list(PCAWGValidID_abs = PCAWGValidID_abs_LCDlist$exposures)
  
  name_list <- names(exposures_list)
  exposures_list <- lapply(seq_along(exposures_list),FUN=function(current_ind){
    current_df <- exposures_list[[current_ind]]
    names(current_df) <- name_list[current_ind]
    return(current_df)
  })
  
  
  names(exposures_list) <- name_list
  
  exposures_list <- lapply(exposures_list, function(list){
    colnames(list) <- colnames(current_id_catalogue_df)
    return(list)
  })
  
  subgroups_df <- data.frame(
    PID = colnames(exposures_list[[1]]), subgroup = "test", sum = 1, 
    compl_sum = 1, index = seq_along(exposures_list[[1]]), col = "#FF0000FF")
  
  complete_df_list <- lapply(name_list, function(current_condition){
    print(current_condition)
    current_exposures_df <- exposures_list[[current_condition]]
    current_sigs <- rownames(current_exposures_df)
    current_signatures_df <- PCAWG_SP_ID_sigs_df[, current_sigs, drop=FALSE]
    suppressWarnings(variateExp(
      in_catalogue_df = current_id_catalogue_df,
      in_sig_df = current_signatures_df,
      in_exposures_df = current_exposures_df,
      in_sigLevel = 0.025, in_delta = 0.4))
  })
  
  complete_df_list <- lapply(complete_df_list, function(current_complete_df){
    current_complete_df$norm_exposure <- current_complete_df$exposure /
      current_complete_df[which(current_complete_df$sig == "total"), "exposure"]
    current_complete_df$norm_lower <- current_complete_df$norm_exposure *
      current_complete_df$relLower
    current_complete_df$norm_upper <- current_complete_df$norm_exposure * 
      current_complete_df$relUpper
    return(current_complete_df)
  })
  
  names(complete_df_list)<-names(exposures_list)
  
  plot_ID <- plotExposuresConfidence_indel(complete_df_list$PCAWGValidID_abs,
                                           subgroups_df,
                                           PCAWG_SP_ID_sigInd_df)
  
  
  return(list(p_complete_PCAWG_ID = plot_ID,
              complete_PCAWG_ID=complete_df_list$PCAWGValidID_abs))
}


#' Plot exposures including confidence intervals for exposures of SNVs and
#' INDELs
#'
#' Plot the exposures to extracted signatures including the confidence intervals
#' computed e.g. by \code{\link[YAPSA]{variateExp}}
#'
#' @param in_complete_df Melted numeric input data frame e.g. as computed by
#'   \code{\link{variateExp}}
#' @param in_subgroups_df Data frame containing meta information on subgroup
#'   attribution of the samples in the cohort of interest.
#' @param in_sigInd_df Data frame with meta information on the signatures used
#'   in the analysis.
#'
#' @return The function returns a gtable object which can be plotted with
#'   \code{plot} or \code{grid.draw}
#'
#' @export
#'
#' @examples NULL
#' 
plotExposuresConfidence_indel <- function(in_complete_df,
                                          in_subgroups_df,
                                          in_sigInd_df) {
  sig_colour_vector <- c("black", in_sigInd_df$colour)
  names(sig_colour_vector) <- c("total", as.character(in_sigInd_df$sig))
  in_complete_df$sample <- factor(in_complete_df$sample, 
                                  levels = in_subgroups_df$PID[order(in_subgroups_df$index)])
  in_complete_df$subgroup <-
    in_subgroups_df$subgroup[match(in_complete_df$sample, 
                                   in_subgroups_df$PID)]
  exposure_plot <- in_complete_df[which(in_complete_df$sig != 
                                          "total"), ] %>% 
    ggplot(aes(x = reorder(sample, -exposure), 
               y = exposure, fill = sig)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
    facet_grid(sig ~ subgroup, scales = "free_x", 
               space = "free_x", switch = "x") + 
    theme_grey() + 
    theme(axis.text.x = element_text(angle = 90,  vjust = 0.5),
          panel.border = element_rect(fill = NA, 
                                      colour = "black"), 
          strip.background = element_rect(
            colour = "black"), 
          legend.position = "none") + 
    scale_fill_manual(values = sig_colour_vector) + 
    scale_x_discrete(name ="sample identifier")
  subgroup_aggregate_df <- aggregate(col ~ subgroup, data = in_subgroups_df, 
                                     FUN = head, 1)
  index_aggregate_df <- aggregate(index ~ subgroup, data = in_subgroups_df, 
                                  FUN = mean)
  subgroup_colour_vector <- 
    subgroup_aggregate_df$col[order(index_aggregate_df$index)]
  exposure_g <- ggplot_gtable(ggplot_build(exposure_plot))
  stripr <- which(grepl("strip-b", exposure_g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl("rect", exposure_g$grobs[[i]]$grobs[[1]]$childrenOrder))
    exposure_g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <-
      subgroup_colour_vector[k]
    k <- k + 1
  }
  total_p <- in_complete_df[which(in_complete_df$sig == "total"),] %>% 
    ggplot(aes(x = reorder(sample, -exposure), y = exposure, fill = sig)) + 
    geom_bar(stat = "identity", fill = "black") + 
    facet_grid(sig ~ subgroup, 
               scales = "free_x",
               space = "free_x", 
               switch = "x") + 
    theme_grey() + theme(axis.text.x = element_blank(), 
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(), 
                         axis.title.y = element_blank(),
                         panel.border = element_rect(fill = NA, 
                                                     colour = "black"), 
                         strip.background = element_rect(colour = "black"), 
                         strip.text.x = element_blank(), 
                         legend.position = "none")
  total_g <- ggplotGrob(total_p)
  all_g <- rbind(total_g, exposure_g)
  return(all_g)
}
