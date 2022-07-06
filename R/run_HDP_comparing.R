#!/usr/bin/env Rscript

#################################################################################################################
####         Expectation Maximisation (EM) to divvy up hdp signatures into known cosimic signatures          #### 
#################################################################################################################
#--> script adapted from https://github.com/HLee-Six/colon_microbiopsies/blob/master/signature_extraction/subsitutions_hdp_signature_extraction/SBS_EM_to_share_out_sigs.R

######## Info ########
# i) For known signatures (e.g. when HDP calls something as “SBS 1”)
#   a.      we calculate the cosine similarity between the HDP version of the signature, and the PCAWG version of the signature.
#   b.      If the cosine similarity is below a certain threshold (0.95 seems to work well), we then apply an expectation maximisation (EM) algorithm to deconvolute this signature into its constituents. E.g. the preconditioned HDP SBS1 goes to 38.5% SBS1, 47.5% SBS5, and 15.0% SBS18. Reassuringly, 38% of the mutations that make up the HDP SBS1 are C to T at CpG, so that fits well.

# ii) For novel signatures (i.e. ones that Nicola’s algorithm calls as novel e.g. "SBS N1")
#   a.      We always apply EM to break an HDP signature down into a composite of PCAWG signatures.
#   b.      We then reconstitute the signature by adding up the PCAWG signatures in the proportion that EM gives us (e.g. IDN1 is 27% ID1, 14% ID2, 20% ID5 etc). I have tried only using signatures that contribute >10% of the mutations to avoid overfitting.
#   c.      We compute the cosine similarity of the reconstituted signature to the known signature. If that is less than a certain threshold, then we consider that the signature is truly novel. Otherwise, we break it down into its constituents and present the data in that way.
######################


#options and libraries
options(stringsAsFactors = F)
set.seed(123)
library(lattice)
library(deconstructSigs)
library(lsa)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(dplyr)
library(cowplot)


#parameters
cmdArgs      <- commandArgs(trailingOnly = TRUE)
hdp_dir      <- cmdArgs[1]
cosineSim_threshold <- as.numeric(cmdArgs[2])
maxiter_EM          <- as.numeric(cmdArgs[3])
EMfrac_threshold    <- as.numeric(cmdArgs[4])
sigRef_file         <- as.character(cmdArgs[5])
cancerSigs_file     <- as.character(cmdArgs[6])
normSigs            <- as.character(cmdArgs[7])


# read in the signature extraction results
if(normSigs == 'TRUE'){
  output_dir  <- paste0(hdp_dir, '/normSignatures/')
  hsbs        <- read.table(paste0(output_dir, '/normSignatures.txt'), sep="\t", header=T)
} else {
  output_dir  <- hdp_dir
  hsbs  <- read.table(paste0(output_dir, '/signatures.txt'), sep="\t", header=T)
}

# read in the pcawg / cosmic sigs
psbs           <- read.table(sigRef_file, header = T)
rownames(psbs) <- psbs[,1]
psbs           <- psbs[match(rownames(hsbs), rownames(psbs)),-1]



#########################################
########        Functions        ########
#########################################

EM_algorithm <- function(output_signature, input_signatures, maxiter = 1000){
  
  num_signatures <- nrow(input_signatures)
  alpha <- runif(num_signatures)
  alpha <- alpha/sum(alpha) # Random start (seems to give ~identical results)
  #alpha = rep(1/num_signatures,num_signatures) # Uniform start
  
  for (iter in 1:maxiter) {
    contr <- t(array(alpha,dim=c(num_signatures,ncol(input_signatures)))) * t(input_signatures)
    probs <- contr/array(rowSums(contr),dim=dim(contr))
    probs[is.na(probs)] <- 0
    probs <- probs * as.numeric(output_signature)
    old_alpha <- alpha
    alpha     <- colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha)) < 1e-5) {
      break
    }
  }
  return(alpha)
}


#function to plot SBS
plot_SBS <- function(signatures_df, title = ''){
  
  #add groups
  colnames(signatures_df) <- c('channel', 'value')
  signatures_df$group <- substr(as.character(signatures_df$channel), start = 3, stop = 5)
  signatures_df <- signatures_df %>%
    mutate(channel = factor(channel, levels = signatures_df$channel),
           group = factor(group, levels = unique(signatures_df$group)),
           value = value * 100)
  
  #set colours
  colours   <- setNames(c('#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'), unique(signatures_df$group))
  strip_name_colours <- c('black','white','white','black','black','black')
  xlabels   <- paste0(substr(as.character(signatures_df$channel), start = 1, stop = 1),
                      substr(as.character(signatures_df$channel), start = 3, stop = 3),
                      substr(as.character(signatures_df$channel), start = 7, stop = 7))
  
  #plot
  p <- ggplot(signatures_df, aes(x = channel, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(.~ group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours, guide = 'none') +
    xlab('') +
    ylab('% SBS') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    scale_x_discrete(labels = xlabels) +
    ggtitle(title) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                            strip.text = element_text(face = 'bold'), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  #change colours and labels of strips
  g       <- ggplot_gtable(ggplot_build(p))
  striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
  k <- 1
  for (i in striprt) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
    t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    
    k <- k+1
  }

  return(g)
}



#function to plot indels
plot_DBS <- function(signatures_df, title = ''){
  
  #add groups
  colnames(signatures_df) <- c('channel', 'value')
  signatures_df$group     <- paste0(matrix(unlist(strsplit(as.character(signatures_df$channel), '>')), ncol = 2, byrow = T)[,1], '>NN')
  signatures_df <- signatures_df %>%
    mutate(channel = factor(channel, levels = signatures_df$channel),
           group = factor(group, levels = unique(signatures_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
                        '#e3211d', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'), unique(signatures_df$group))
  strip_name_colours <- c('black','white','black','white','black','white','black','white','black','white')
  xlabels <- matrix(unlist(strsplit(as.character(signatures_df$channel), '>')), ncol = 2, byrow = T)[,2]
  
  #plot
  p <- ggplot(signatures_df, aes(x = channel, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(.~ group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours, guide = 'none') +
    xlab('') +
    ylab('% DBS') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    scale_x_discrete(labels = xlabels) +
    ggtitle(title) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                            strip.text = element_text(face = 'bold'), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  #change colours and labels of strips
  g       <- ggplot_gtable(ggplot_build(p))
  striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
  k <- 1
  for (i in striprt) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
    t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    
    k <- k+1
  }
  
  return(g)
}



#function to plot indels
plot_ID <- function(signatures_df, title = ''){
  
  #add groups
  colnames(signatures_df) <- c('channel', 'value')
  signatures_df$group <- paste(matrix(unlist(strsplit(as.character(signatures_df$channel), ':')), ncol = 4, byrow = T)[,1],
                               matrix(unlist(strsplit(as.character(signatures_df$channel), ':')), ncol = 4, byrow = T)[,2],
                               matrix(unlist(strsplit(as.character(signatures_df$channel), ':')), ncol = 4, byrow = T)[,3], sep = ':')
  signatures_df <- signatures_df %>%
    mutate(channel = factor(channel, levels = signatures_df$channel),
           group = factor(group, levels = unique(signatures_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#fdbe6f', '#ff8002', '#b0dd8b', '#36a12e', '#fdcab5', '#fc8a6a', '#f14432','#bc191a', '#d0e1f2', '#94c4df', '#4a98c9', '#1764ab', 
                        '#e1e1ef', '#b6b6d8', '#8683bd','#62409b'), unique(signatures_df$group))
  legend_labels <- setNames(c('1bp C Deletion', '1bp T Deletion','1bp C Insertion', '1bp T Insertion',
                              '2bp Deletion at Repeats', '3bp Deletion at Repeats', '4bp Deletion at Repeats', '5+bp Deletion at Repeats',
                              '2bp Insertion at Repeats', '3bp Insertion at Repeats', '4bp Insertion at Repeats', '5+bp Insertion at Repeats',
                              '2bp Microhomology Deletion', '3bp Microhomology Deletion', '4bp Microhomology Deletion', '5+bp Microhomology Deletion'),
                            unique(signatures_df$group))
  strip_name_colours <- rep('black', length(colours))
  strip_name_colours[which(names(colours) %in% c("1:Del:T", "1:Ins:T", "5:Del:R", "5:Ins:R", "5:Del:M"))] <- 'white'
  strip_name <- c('C', 'T', 'C', 'T', '2', '3', '4', '5+', '2', '3', '4', '5+', '2', '3', '4', '5+')
  
  xlabels <- c('1','2','3','4','5', '6+', '1','2','3','4','5', '6+',
               '0','1','2','3','4','5+','0','1','2','3','4','5+',
               '1','2','3','4','5', '6+', '1','2','3','4','5', '6+',
               '1','2','3','4','5', '6+', '1','2','3','4','5', '6+',
               '0','1','2','3','4','5+','0','1','2','3','4','5+',
               '0','1','2','3','4','5+','0','1','2','3','4','5+',
               '1','1','2','1','2','3','1','2','3','4','5+')
  
  #plot
  p <- ggplot(signatures_df, aes(x = channel, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(.~group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours, labels = legend_labels) +
    xlab('Homopolymer Length / Number of Repeat Units / Microhomology Length') +
    ylab('% Indels') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    scale_x_discrete(breaks = as.character(signatures_df$channel), labels = xlabels) +
    ggtitle(title) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15),
                            legend.position = 'none', strip.text = element_text(face = 'bold'),
                            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  #change colours and labels of strips
  g       <- ggplot_gtable(ggplot_build(p))
  striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
  k <- 1
  for (i in striprt) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
    t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$label  <- strip_name[k]
    
    k <- k+1
  }
  
  return(g)
}


#function to plot SVs
plot_SV <- function(signatures_df, title = ''){
  
  #add groups
  colnames(signatures_df) <- c('channel', 'value')
  signatures_df$group <- sapply(as.character(signatures_df$channel), function(x){ unlist(strsplit(x, ':'))[1] })
  signatures_df <- signatures_df %>%
    mutate(channel = factor(channel, levels = signatures_df$channel),
           group = factor(group, levels = unique(signatures_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#1f78b4', '#33a02c', '#e31a1c','#ff7f00'), unique(signatures_df$group))
  strip_name_colours <- c(rep('black', 4), rep('white', 4))
  strip_name <- matrix(unlist(strsplit(as.character(unique(signatures_df$group)), '_')), ncol = 2, byrow = T)[,2]
  
  legend_labels <- setNames(c('Clustered Deletion', 'Clustered Tandem Duplication','Clustered Inversion', 'Clustered Translocation',
                              'Non-clustered Deletion', 'Non-clustered Tandem Duplication','Non-clustered Inversion', 'Non-clustered Translocation'),
                            unique(signatures_df$group))
  
  xlabels <- as.character(sapply(as.character(signatures_df$channel), function(x){ unlist(strsplit(x, ':'))[2] }))
  xlabels <- xlabels[!is.na(xlabels)] 

  #plot
  p <- ggplot(signatures_df, aes(x = channel, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(.~group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours[c(1,5,2,6,3,7,4,8)], labels = legend_labels[c(1,5,2,6,3,7,4,8)]) +
    xlab('Rearrangment size') +
    ylab('% SVs') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    scale_x_discrete(breaks = grep('trans', as.character(signatures_df$channel), invert = T, value = T), labels = xlabels) +
    ggtitle(title) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                            strip.text = element_text(face = 'bold'), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), legend.position = 'none')
  
  #change colours and labels of strips
  g       <- ggplot_gtable(ggplot_build(p))
  striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
  k <- 1
  for (i in striprt) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
    t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$label  <- strip_name[k]
    
    k <- k+1
  }
  
  return(g)
}


#function to plot SCNAs
plot_SCNA <- function(signatures_df, title = ''){
  
  #add groups
  colnames(signatures_df) <- c('channel', 'value')
  signatures_df$CN        <- lapply(as.character(signatures_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[1] })
  signatures_df$type      <- lapply(as.character(signatures_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[2] })
  signatures_df$length    <- lapply(as.character(signatures_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[3] })
  signatures_df$value     <- signatures_df$value * 100
  signatures_df$group     <- paste0(signatures_df$CN, ":", signatures_df$type)
  signatures_df$group     <- factor(signatures_df$group, levels = unique(signatures_df$group))
  
  signatures_df$type      <- factor(signatures_df$type, levels = c("homdel", "LOH", "het"))
  signatures_df$length    <- as.character(signatures_df$length)
  signatures_df$channel   <- factor(signatures_df$channel, levels = unique(signatures_df$channel))
  
  
  #set colours
  colours       <- setNames(c('#fb8072', '#ffd92f','#66c2a5', '#e78ac3',  '#8da0cb', '#fc8d62', '#1b9e77', '#e7298a',  '#7570b3', '#d95f02'), unique(signatures_df$group))
  legend_labels <- setNames(c('homozygous deletion', 'LOH & CN1', 'LOH & CN2', 'LOH & CN3-4', 'LOH & CN5-8', 'LOH & CN9+', 
                              'heterozygous & CN2', 'heterozygous & CN3-4', 'heterozygous & CN5-8', 'heterozygous & CN9+'), unique(signatures_df$group))
  
  strip_name_colours <- rep('black', length(colours))
  strip_name_colours[which(names(colours) %in% c("2:het", "3-4:het", "5-8:het", "9+:het"))] <- 'white'
  strip_name <- c('CN 0', 'CN 1', 'CN 2', 'CN 3-4', 'CN 5-8', 'CN 9+', 'CN 2', 'CN 3-4', 'CN 5-8', 'CN 9+')
  
  xlabels <- c("0-100kb", "100kb-1Mb", ">1Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb",
               "0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb")
  #plot
  p <- ggplot(signatures_df, aes(x = channel, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(.~group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours[c(1,2,3,7,4,8,5,9,6,10)], labels = legend_labels[c(1,2,3,7,4,8,5,9,6,10)]) +
    xlab('Segment Length') +
    ylab('% CN segments') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    scale_x_discrete(breaks = signatures_df$channel, labels = xlabels) +
    ggtitle(title) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = 'bottom', strip.text = element_text(face = 'bold'),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  #change colours and labels of strips
  g       <- ggplot_gtable(ggplot_build(p))
  striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
  k <- 1
  for (i in striprt) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
    t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$label  <- strip_name[k]
    
    k <- k+1
  }
  
  return(g)
}





#########################################
########           Main          ########
#########################################

novel <- grep("N", colnames(hsbs), value=T)
known <- grep("N", colnames(hsbs), value=T, invert = T)
# I think leave the residual as it is.

# calculate cosine similarity for all of them
cosmat           <- data.frame(matrix(nrow=ncol(hsbs), ncol=ncol(psbs)))
rownames(cosmat) <- colnames(hsbs)
colnames(cosmat) <- colnames(psbs)

for (i in 1:nrow(cosmat)) {
  for (j in 1:ncol(cosmat)) {
    cosmat[i,j] <- cosine(x=hsbs[,rownames(cosmat)[i]], y=psbs[,colnames(cosmat)[j]])
  }
}

colnames(cosmat) <- paste0("cosmic_", colnames(cosmat))
rownames(cosmat) <- paste0("hdp_", rownames(cosmat))
write.table(cosmat, file = paste0(output_dir, "cosineSimilarities.txt"), sep="\t", col.names = T, row.names = T, quote=F)

#plot
plot_data <- reshape2::melt(as.matrix(cosmat))
colnames(plot_data) <- c('HDP', 'COMSIC', 'cosineSimilarity')
plot_data$HDP    <- sub('hdp_', '', plot_data$HDP)
plot_data$COMSIC <- sub('cosmic_', '', plot_data$COMSIC)
plot_data$COMSIC <- factor(plot_data$COMSIC, levels = unique(plot_data$COMSIC))

p <- ggplot(plot_data, aes(x = COMSIC, y = HDP, fill = cosineSimilarity)) + 
  geom_tile(colour = 'black') +
  scale_fill_gradient(low = '#fee090', high = '#a50026') +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_text(aes(label = ifelse(cosineSimilarity > cosineSim_threshold, '*', '')), colour = 'white') +
  labs(title = 'Cosine Similarities',
       subtitle = paste0('(Cosine Similarity cutoff = ', cosineSim_threshold, ')')) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))

pdf(paste0(output_dir, "cosineSimilarities.pdf"), height=4, width=15)
plot(p)
dev.off()



###############
#### known signatures ####
##############

# only split it into signatures previously found in that cancer type if provided
signatures <- t(psbs)
prefix     <- ''
if(cancerSigs_file != 'NA'){
  cancer_sigs <- readRDS(cancerSigs_file)
  signatures  <- signatures[cancer_sigs,]
  prefix      <- 'cancerSigs_'
}


# for known sigs, look at the cosine similarity between the hdp and the cosmic version and check if signatures need to be split
if(length(known) > 0){
  cosmat[paste0("hdp_", known),paste0("cosmic_", known)]
  check_sigs <- known[diag(as.matrix(cosmat[paste0("hdp_", known),paste0("cosmic_", known)])) < cosineSim_threshold]
} else {
  check_sigs <- c()
}

if(length(check_sigs) > 0){
  mutations <- hsbs[,check_sigs, drop = F]
  
  #run EM-algorithm
  signature_fraction <- lapply(check_sigs, function(x){
    output_signature <- mutations[,x]
    output_signature[is.na(output_signature)] <- 0
    as.matrix(EM_algorithm(output_signature, signatures, maxiter_EM))
  })
  signature_fraction <- Reduce(cbind, signature_fraction)
  colnames(signature_fraction) <- check_sigs
  
  #plot
  pdf(paste0(output_dir, prefix, "EM_breakup_knownSigs.pdf"), height = 4, width = 15)
  color.palette <- colorRampPalette(c("white", "orange", "purple"))
  print(levelplot((signature_fraction[dim(signature_fraction)[1]:1,,drop = F]),
                  col.regions=color.palette, seq(0,1, 0.1),
                  aspect="fill", scales=list(x=list(rot=45), tck = c(1,0)),
                  xlab = 'COSMIC', ylab = 'HDP', main = 'EM breakp-up (%) of known signatures'))
  dev.off()
  
  #save results
  s           <- signature_fraction
  rownames(s) <- paste0('cosmic_',rownames(s))
  colnames(s) <- paste0('hdp_',colnames(s))
  write.table(s, file = paste0(output_dir, prefix, "EM_breakup_knownSigs.txt"), sep="\t", col.names=T, row.names = T, quote=F)
  
  #re-run EM only with signatures with fraction greater than EMfrac_threshold and cosine similarity > 0.8
  subset_signature_fractions <- lapply(check_sigs, function(x){
    print(x)
    constit     <- names(which(signature_fraction[,x] > EMfrac_threshold))
    if(length(sub('cosmic_', '', colnames(cosmat)[which(cosmat[paste0('hdp_', x),] >= cosineSim_threshold)])) == 1){
      out <- data.frame(hdp = x, cosmic = sub('cosmic_', '', colnames(cosmat)[which(cosmat[paste0('hdp_', x),] >= cosineSim_threshold)]), fraction = 1)
      return(out)
    }
    cosine_sigs <- sub('cosmic_', '', colnames(cosmat)[which(cosmat[paste0('hdp_', x),] > 0.8)])
    
    sig_pairs   <- t(combn(unique(c(constit, cosine_sigs)), 2))
    pairwise_test <- lapply(1:nrow(sig_pairs), function(i){
      sigs        <- as.matrix(psbs[,sig_pairs[i,], drop = F])
      input_signatures <- t(sigs)
      output_signature <- hsbs[,x]
      results    <- EM_algorithm(output_signature, input_signatures, maxiter_EM)
      # while(any(results <= EMfrac_threshold)){
      #   constit     <- names(which(results > EMfrac_threshold))
      #   sigs        <- as.matrix(psbs[,constit, drop = F])
      #   input_signatures <- t(sigs)
      #   output_signature <- hsbs[,x]
      #   results    <- EM_algorithm(output_signature, input_signatures, maxiter_EM)
      # }
      
      reconstr <- results[1] * psbs[,sig_pairs[i,1]] + results[2] * psbs[,sig_pairs[i,2]]
      cosSim   <- cosine(reconstr, output_signature)
      out      <- data.frame(hdp = x, cosmic = names(results), fraction = as.numeric(results), cosineSimilarity = cosSim)
      return(out)
    })
    
    cosineSims <- sapply(pairwise_test, function(y) y$cosineSimilarity[1])
    index      <- which(cosineSims == max(cosineSims))
    
    out <- pairwise_test[[index]]
    out <- out[, -1*ncol(out)]
    return(out)
  })
  subset_signature_fractions <- Reduce(rbind, subset_signature_fractions)  
  
  fraction_matrix <- matrix(NA, nrow = length(check_sigs), ncol = length(unique(subset_signature_fractions$cosmic)),
                            dimnames = list(check_sigs, unique(subset_signature_fractions$cosmic)))
  for(x in check_sigs){
    sub <- subset_signature_fractions[subset_signature_fractions$hdp == x,]
    fraction_matrix[x, sub$cosmic] <- sub$fraction
  }
  subset_signature_fractions <- t(fraction_matrix)
  subset_signature_fractions <- subset_signature_fractions[order(as.numeric(sub('SBS', '', rownames(subset_signature_fractions)))),,drop = F]
  
  #plot
  pdf(paste0(output_dir, prefix, "EM_breakup_knownSigs_signif.pdf"), height = 4, width = 6)
  color.palette <- colorRampPalette(c("white", "orange", "purple"))
  print(levelplot((subset_signature_fractions[dim(subset_signature_fractions)[1]:1,,drop = F]),
                  col.regions=color.palette, at = seq(0,1, 0.1),
                  aspect="fill", scales=list(x=list(rot=45), tck = c(1,0)),
                  xlab = 'COSMIC', ylab = 'HDP', main = 'EM breakp-up (%) of known signatures'))
  dev.off()
  
  #save results
  s           <- subset_signature_fractions
  rownames(s) <- paste0('cosmic_',rownames(s))
  colnames(s) <- paste0('hdp_',colnames(s))
  write.table(s, file = paste0(output_dir, prefix, "EM_breakup_knownSigs_signif.txt"), sep="\t", col.names=T, row.names = T, quote=F)
  
  
  
  # reconstitute HDP signatures using the identified COSMIC sigs and calculate the cosine similarity between the original and the reconstituted sig.
  # --> skip if only one COMSIC signature is assigned
  reconSBS <- lapply(check_sigs, function(x){
    constit  <- names(which(!is.na(subset_signature_fractions[,x])))
    if(length(constit) == 1 && constit == x){
      return(NULL)
    }
    reconsbs <- rep(0, 96)
    for (ct in constit) {
      reconsbs <- reconsbs + (psbs[,ct]*subset_signature_fractions[ct,x])
    }
    reconsbs <- data.frame(reconsbs)
    colnames(reconsbs) <- x
    return(t(reconsbs))
  })
  reconSBS <- Reduce(rbind, reconSBS)
  
  cosineSimilarities <- sapply(rownames(reconSBS), function(x) cosine(reconSBS[x,], hsbs[,x]))
  
  ### Plot ###
  # 1) plot the original
  # 2) their broken down signatures
  # 3) the reconstituted sigs
  
  pdf(paste0(output_dir, prefix, "HDP_reconstitution_knownSigs.pdf"), width = ifelse(nrow(psbs) < 32, 9, 11), height = 12)
  for(i in rownames(reconSBS)){
    constit <- names(which(!is.na(subset_signature_fractions[,i])))
    plot_sigs <- data.frame(channel = rownames(psbs), original = hsbs[,i], reconstructed = reconSBS[i,])
    for (ct in constit) {
      plot_sigs <- cbind(plot_sigs, psbs[,ct])
      colnames(plot_sigs)[ncol(plot_sigs)] <- paste0(ct, " accounts for ", round(subset_signature_fractions[ct,i], digits=2))
    }
    y         <- c(0.75, 0.5, 0.25,0)
    full_plot <- ggdraw()
    for(x in 2:ncol(plot_sigs)){
      title <- ifelse(colnames(plot_sigs)[x] == 'original', paste0('Original ', i), colnames(plot_sigs)[x])
      title <- ifelse(title == 'reconstructed', paste0('Reconstituted ', i, ", cosine similarity to original: ", round(cosineSimilarities[i], digits=2)), title)
      if(nrow(psbs) == 96){ p <- plot_SBS(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 78){ p <- plot_DBS(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 83){ p <- plot_ID(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 32){ p <- plot_SV(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 48){ p <- plot_SCNA(plot_sigs[,c(1,x)], title) }
      full_plot <- full_plot + draw_plot(p, x = 0, y = y[x-1], width = 1, height = 0.25)
    }
    plot(full_plot)
  }
  dev.off()
}



####################
######## novel signatures
####################
# now break down all novel signatures, rebuild them back up, and see how well each is explained.
# only split it into the good signatures.
check_sigs <- novel

if(length(check_sigs) > 0){
  mutations <- hsbs[,check_sigs, drop = F]
  
  #run EM-algorithm
  signature_fraction <- lapply(check_sigs, function(x){
    output_signature <- mutations[,x]
    output_signature[is.na(output_signature)] <- 0
    as.matrix(EM_algorithm(output_signature, signatures, maxiter_EM))
  })
  signature_fraction <- Reduce(cbind, signature_fraction)
  colnames(signature_fraction) <- check_sigs
  
  #plot
  pdf(paste0(output_dir, prefix, "EM_breakup_novelSigs.pdf"), height = 4, width = 15)
  color.palette <- colorRampPalette(c("white", "orange", "purple"))
  print(levelplot((signature_fraction[dim(signature_fraction)[1]:1,,drop = F]),
                  col.regions=color.palette, seq(0,1, 0.1),
                  aspect="fill", scales=list(x=list(rot=45), tck = c(1,0)),
                  xlab = 'COSMIC', ylab = 'HDP', main = 'EM breakp-up (%) of novel signatures'))
  dev.off()
  
  #save results
  s           <- signature_fraction
  rownames(s) <- paste0('cosmic_',rownames(s))
  colnames(s) <- paste0('hdp_',colnames(s))
  write.table(s, file = paste0(output_dir, prefix, "EM_breakup_novelSigs.txt"), sep="\t", col.names=T, row.names = T, quote=F)
  
  
  #re-run EM only with signatures with fraction greater than EMfrac_threshold and cosine similarity > 0.8
  subset_signature_fractions <- lapply(check_sigs, function(x){
    print(x)
    constit     <- names(which(signature_fraction[,x] > EMfrac_threshold))
    if(length(sub('cosmic_', '', colnames(cosmat)[which(cosmat[paste0('hdp_', x),] >= cosineSim_threshold)])) > 0){
      if(length(sub('cosmic_', '', colnames(cosmat)[which(cosmat[paste0('hdp_', x),] >= cosineSim_threshold)])) > 1){
        cosmic <- cosmat[paste0('hdp_', x),which(cosmat[paste0('hdp_', x),] >= cosineSim_threshold)]
        cosmic <- sub('cosmic_', '', colnames(cosmic)[which.max(cosmic)])
      } else {
        cosmic <- sub('cosmic_', '', colnames(cosmat)[which(cosmat[paste0('hdp_', x),] >= cosineSim_threshold)])
      }
      out <- data.frame(hdp = x, cosmic = cosmic, fraction = 1)
      return(out)
    }
    cosine_sigs <- sub('cosmic_', '', colnames(cosmat)[which(cosmat[paste0('hdp_', x),] > 0.8)])
    
    if(length(constit) <= 1){
      constit <- names(sort(signature_fraction[,x], decreasing = T)[1:2])
    }
    
    sig_pairs     <- t(combn(unique(c(constit, cosine_sigs)), 2))
    pairwise_test <- lapply(1:nrow(sig_pairs), function(i){
      sigs        <- as.matrix(psbs[,sig_pairs[i,], drop = F])
      input_signatures <- t(sigs)
      output_signature <- hsbs[,x]
      results    <- EM_algorithm(output_signature, input_signatures, maxiter_EM)
      # while(any(results <= EMfrac_threshold)){
      #   constit     <- names(which(results > EMfrac_threshold))
      #   sigs        <- as.matrix(psbs[,constit, drop = F])
      #   input_signatures <- t(sigs)
      #   output_signature <- hsbs[,x]
      #   results    <- EM_algorithm(output_signature, input_signatures, maxiter_EM)
      # }
      
      reconstr <- results[1] * psbs[,sig_pairs[i,1]] + results[2] * psbs[,sig_pairs[i,2]]
      cosSim   <- cosine(reconstr, output_signature)
      out      <- data.frame(hdp = x, cosmic = names(results), fraction = as.numeric(results), cosineSimilarity = cosSim)
      return(out)
    })
    
    cosineSims <- sapply(pairwise_test, function(y) y$cosineSimilarity[1])
    index      <- which(cosineSims == max(cosineSims))
    
    out <- pairwise_test[[index]]
    out <- out[, -1*ncol(out)]
    return(out)
  })
  subset_signature_fractions <- Reduce(rbind, subset_signature_fractions)  
  
  fraction_matrix <- matrix(NA, nrow = length(check_sigs), ncol = length(unique(subset_signature_fractions$cosmic)),
                            dimnames = list(check_sigs, unique(subset_signature_fractions$cosmic)))
  for(x in check_sigs){
    sub <- subset_signature_fractions[subset_signature_fractions$hdp == x,]
    fraction_matrix[x, sub$cosmic] <- sub$fraction
  }
  subset_signature_fractions <- t(fraction_matrix)
  subset_signature_fractions <- subset_signature_fractions[order(as.numeric(sub('SBS|DBS|ID', '', rownames(subset_signature_fractions)))),,drop = F]
  
  #plot
  pdf(paste0(output_dir, prefix, "EM_breakup_novelSigs_signif.pdf"), height = 4, width = 6)
  color.palette <- colorRampPalette(c("white", "orange", "purple"))
  print(levelplot((subset_signature_fractions[dim(subset_signature_fractions)[1]:1,,drop = F]),
                  col.regions=color.palette, at = seq(0,1, 0.1),
                  aspect="fill", scales=list(x=list(rot=45), tck = c(1,0)),
                  xlab = 'COSMIC', ylab = 'HDP', main = 'EM breakp-up (%) of known signatures'))
  dev.off()
  
  #save results
  s           <- subset_signature_fractions
  rownames(s) <- paste0('cosmic_',rownames(s))
  colnames(s) <- paste0('hdp_',colnames(s))
  write.table(s, file = paste0(output_dir, prefix, "EM_breakup_novelSigs_signif.txt"), sep="\t", col.names=T, row.names = T, quote=F)
  
  
  
  # reconstitute HDP signatures using the identified COSMIC sigs and calculate the cosine similarity between the original and the reconstituted sig.
  # --> skip if only one COMSIC signature is assigned
  reconSBS <- lapply(check_sigs, function(x){
    constit  <- names(which(!is.na(subset_signature_fractions[,x])))
    if(length(constit) == 1 && constit == x){
      return(NULL)
    }
    reconsbs <- rep(0, nrow(mutations))
    for (ct in constit) {
      reconsbs <- reconsbs + (psbs[,ct]*subset_signature_fractions[ct,x])
    }
    reconsbs <- data.frame(reconsbs)
    colnames(reconsbs) <- x
    return(t(reconsbs))
  })
  reconSBS <- Reduce(rbind, reconSBS)
  
  cosineSimilarities <- sapply(rownames(reconSBS), function(x) cosine(reconSBS[x,], hsbs[,x]))
  
  ### Plot ###
  # 1) plot the original
  # 2) their broken down signatures
  # 3) the reconstituted sigs
  
  pdf(paste0(output_dir, prefix, "HDP_reconstitution_novelSigs.pdf"), width = ifelse(nrow(psbs) < 50, 9, 11), height = 12)
  for(i in rownames(reconSBS)){
    constit <- names(which(!is.na(subset_signature_fractions[,i])))
    plot_sigs <- data.frame(channel = rownames(psbs), original = hsbs[,i], reconstructed = reconSBS[i,])
    for (ct in constit) {
      plot_sigs <- cbind(plot_sigs, psbs[,ct])
      colnames(plot_sigs)[ncol(plot_sigs)] <- paste0(ct, " accounts for ", round(subset_signature_fractions[ct,i], digits=2))
    }
    y         <- c(0.75, 0.5, 0.25,0)
    full_plot <- ggdraw()
    for(x in 2:ncol(plot_sigs)){
      title <- ifelse(colnames(plot_sigs)[x] == 'original', paste0('Original ', i), colnames(plot_sigs)[x])
      title <- ifelse(title == 'reconstructed', paste0('Reconstituted ', i, ", cosine similarity to original: ", round(cosineSimilarities[i], digits=2)), title)
      if(nrow(psbs) == 96){ p <- plot_SBS(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 78){ p <- plot_DBS(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 83){ p <- plot_ID(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 32){ p <- plot_SV(plot_sigs[,c(1,x)], title) }
      if(nrow(psbs) == 48){ p <- plot_SCNA(plot_sigs[,c(1,x)], title) }
      full_plot <- full_plot + draw_plot(p, x = 0, y = y[x-1], width = 1, height = 0.25)
    }
    plot(full_plot)
  }
  dev.off()
}  








