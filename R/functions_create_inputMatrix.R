##################################################################################################
########            Functions to create different types of input matrices for HDP         ########
##################################################################################################


### SBS signature input matrix ###
#--> use deconstructSigs functions 
create_SBS_inputMatrix <- function(mutTable, bsg = NULL){
  #first 4 columns need to be in order of Sample_id, chr, pos, ref and alt
  input <- mutTable[,1:5]
  colnames(input) <- c('Sample', 'chr', 'pos', 'ref', 'alt')
  input$pos       <- as.numeric(as.character(input$pos))
  input_matrix <- deconstructSigs::mut.to.sigs.input(input, bsg = bsg)
  return(input_matrix)
}


### DBS signature input matrix ###
#--> use deconstructSigs functions 
create_DBS_inputMatrix <- function(mutTable){
  #first 5 columns need to be in order of Sample_id, chr, pos, ref, alt and alt_count
  input <- mutTable[,1:6]
  colnames(input) <- c('Sample', 'chr', 'pos', 'ref', 'alt', 'alt_count')
  input$pos       <- as.numeric(as.character(input$pos))
  
  #check if DBS have been classified in the data already
  if(!any(nchar(input$ref) == 2 & nchar(input$alt) == 2)) {
    DBS_table <- lapply(unique(input$Sample), function(x){
      sample_table <- input[input$Sample == x,]
      sample_table <- sample_table[order(sample_table$chr, sample_table$pos, decreasing = T),]
      ww <- which(diff(sample_table$pos) == 1)
      if(length(ww) == 0){
        return(NULL)
      } 
      DBS <- lapply(ww, function(i){
        data <- sample_table[c(i, i+1),] 
        if(data$alt_count[1] != data$alt_count[2]){
          return(NULL)
        }
        data.frame(Sample = x, chr = data$chr[1], pos = data$pos[1],ref = paste(data$ref, collapse = ''), alt = paste(data$alt, collapse = ''))
      })
      DBS <- Reduce(rbind, DBS)  
      return(DBS)
    })
    input <- Reduce(rbind, DBS_table)
  }
    
  #get DBS input matrix
  input_matrix <- deconstructSigs::mut.to.sigs.input(input, sig.type = 'DBS')
  return(input_matrix)
}
  
  


### INDEL signature input matrix (83 channels) ###
create_INDEL_inputMatrix <- function(mutTable){
  
  source('/camp/lab/swantonc/working/dietzem/Functions/HDP/R/YAPSA_indel_functions.R')
  
  #first 4 columns need to be in order of Sample_id, chr, pos, ref and alt
  input <- mutTable[,1:5]
  colnames(input) <- c('Sample', 'chr', 'pos',  'ref', 'alt')
  input$pos       <- as.numeric(as.character(input$pos))
  
  #fix indels if needed
  if(any(grep("-", input$alt))){
    input[grep("-", input$alt), "ref"] <- paste0(input[grep("-", input$alt), "ref"], sub('-', '', input[grep("-", input$alt), "alt"]))
    input[grep("-", input$alt), "alt"] <- substr(input[grep("-", input$alt), "ref"],1,1)
  }
  if(any(grep("-", input$ref))){
    ww <- grep("-", input$ref)
    input[ww, "pos"] <- input[ww, "pos"] - 1
    insertion_gr <- GRanges(seqnames = input$chr[ww], IRanges(start = input$pos[ww], end = input$pos[ww]))
    ref_seq      <- getSeq(BSgenome.Hsapiens.UCSC.hg19,insertion_gr)
    input[ww, "ref"] <- as.data.frame(ref_seq)[,1]
    input[ww, "alt"] <- paste0(input[ww, "ref"], input[ww, "alt"])
  }
  
  #filter for indels
  input <- input[nchar(input$ref) > 1 | nchar(input$alt) > 1,]
  if(any(nchar(input$ref) == 2 & nchar(input$alt) == 2)){
    input <- input[!(nchar(input$ref) == 2 & nchar(input$alt) == 2),]
  }
  
  #create vcf like table
  vcf_input <- input[, c('chr', 'pos',  'ref', 'alt', 'Sample')]
  colnames(vcf_input) <- c('CHROM', 'POS', 'REF', 'ALT', 'PID')
  
  #get indel input matrix
  input_matrix <- create_indel_mutation_catalogue_from_df(vcf_input)
  return(input_matrix)
  
}



### SV signature input matrix ###
#--> based on method from Nik-Zainail (https://github.com/Nik-Zainal-Group/signature.tools.lib)












### SCNA signature input matrix ###
#--> based on method from Chris Steel 2021
copyNumberInfo <-  function(chrom, start, end, CN, CNb = NULL, doVariables = TRUE, returnSep = FALSE){
  #segment lengths
  print("lengths")
  lengths = (end - start) / 1000000
  
  #allelic imbalance
  print("imbalance")
  CNa = CN - CNb
  
  #LOH
  LOH = pmin(CNa, CNb)
  
  #combine
  print("combine")
  if(!doVariables){
    combined = list(CN = CN, lengths = lengths, LOH = LOH)
    return(combined)
  }
  
  LOHstatus = ifelse(LOH == 0, "LOH", "het")
  LOHstatus[which(CN == 0)] = "homdel"
  
  varCN = cut(CN,
              breaks=c(-0.5, 0.5, 1.5, 2.5, 4.5, 8.5, Inf),
              labels=c("0", "1", "2","3-4","5-8","9+"))
  
  varLength = rep(NA,length=length(varCN))
  
  hdIndex = which(LOHstatus=="homdel")
  
  if(length(hdIndex) > 0){
    varLength[hdIndex] = paste0(cut(lengths[hdIndex], breaks=c(-0.01, 0.1, 1, Inf)))
    varLength[-hdIndex] = paste0(cut(lengths[-hdIndex], breaks=c(-0.01, 0.1, 1, 10, 40, Inf)))
  } else {
    varLength = paste0(cut(lengths,breaks=c(-0.01, 0.1, 1, 10, 40, Inf)))
  }
  
  renameVarLengths = c("(-0.01, 0.1]" = "0-100kb", "(0.1,1]" = "100kb-1Mb", "(1,10]" = "1Mb-10Mb", 
                       "(10,40]" = "10Mb-40Mb","(40,Inf]" = ">40Mb", "(1,Inf]" = ">1Mb")
  
  varLength = renameVarLengths[paste0(varLength)]
  
  sepVars = paste(varCN, LOHstatus, varLength, sep = ":")
  
  if(returnSep){ return(sepVars) }
  
  variables = table(sepVars)
  return(variables)
}







