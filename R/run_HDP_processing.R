#!/usr/bin/env Rscript

###################################################################################
######                Combine iterations and process results                 ######
###################################################################################
#--> script adapted from Oriol and is partially based on Nicola Robert's examples

# options and libraries
options(stringsAsFactors = F)
library(hdp)
library(RColorBrewer)
library(ggplot2)
library(lsa)
library(Metrics)
library(dplyr)



#parameters
cmdArgs <- commandArgs(trailingOnly = TRUE)
input_matrix_file     <- cmdArgs[1]
output_dir            <- cmdArgs[2]
treeLayer_file        <- cmdArgs[3]
priors_file           <- cmdArgs[4]             #NA = no priors
nMut_cutoff           <- as.numeric(cmdArgs[5]) #NA = no minimum mutations cutoff
sigActivity_cutoff    <- as.numeric(cmdArgs[6])
cohort_cutoff         <- as.numeric(cmdArgs[7])


iteration_dir <- paste0(output_dir, '/iterations/')


#####################################
#####        Functions          #####
#####################################

#plot exposures as boxplot with activity and cohort thresholds
plot_hdp_exposure_boxplot <- function(hdpsample, dpindices, component_names = NULL, sig_active_cutoff = 0.1, cohort_threshold = 0.05){
  
  # input checks
  if (!class(hdpsample) %in% c("hdpSampleChain", "hdpSampleMulti")) {
    stop("hdpsample must have class hdpSampleChain or hdpSampleMulti")
  }
  # if (!validObject(hdpsample)) stop("hdpsample not valid") # too slow on big objects
  if (length(comp_categ_counts(hdpsample)) == 0) {
    stop("No component info for hdpsample. First run hdp_extract_components")
  }
  
  dp_distn <- comp_dp_distn(hdpsample)
  ndp      <- nrow(dp_distn$mean)
  ncomp    <- ncol(dp_distn$mean)
  
  # mean exposures
  exposures <- t(dp_distn$mean[dpindices,,drop=FALSE])
  
  # exclude non-significant exposures
  cis <- dp_distn$cred.int[dpindices]
  nonsig <- lapply(cis, function(x) which(x[1,]==0))
  for (i in 1:length(nonsig)){
    exposures[nonsig[[i]],i] <- 0
  }
  
  #add component_names if provided
  if(!is.null(component_names)){
    rownames(exposures) <- component_names
  }
  
  # exclude component 0
  exposures <- exposures[-1,]
  
  # prepare data for plotting
  plot_data           <- reshape2::melt(exposures)
  colnames(plot_data)[1:2] <- c('Component', 'Sample')
  plot_data$sig_active <- plot_data$value > sig_active_cutoff
  
  cohort_threshold_number <- cohort_threshold * ncol(exposures)
  sig_cohort <- lapply(unique(plot_data$Component), function(x){
    sub            <- plot_data[plot_data$Component == x,]
    sub$sig_cohort <- sum(sub$sig_active) > cohort_threshold_number
    return(sub)
  })
  plot_data <- Reduce(rbind, sig_cohort)
  
  order_signatures <- sapply(unique(plot_data$Component), function(x){
    sub <- plot_data[plot_data$Component == x,]
    sum(sub$sig_active)
  })
  plot_data$Component <- factor(plot_data$Component, levels = unique(plot_data$Component)[order(order_signatures, decreasing = T)])
  
  # plot
  p <- ggplot(plot_data, aes(x = Component, y = value)) + 
    geom_jitter(aes(colour = sig_active, shape = sig_cohort), width = 0.2) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
    scale_colour_manual(name = 'sigActivity PASS', values = c('FALSE' = '#f4a582', 'TRUE' = '#b2182b')) +
    scale_shape_manual(name = 'cohort PASS', values = c('FALSE' = 4, 'TRUE' = 16)) +
    xlab('Component') + ylab('% Activity') + 
    labs(title = 'Signature Activity',
         subtitle = paste0('(sigActivity cutoff = ', sig_active_cutoff, ', cohort cutoff = ', cohort_threshold, ')')) +
    theme_bw()
  
  #fail signatures
  exclude_components <- plot_data[!duplicated(plot_data$Component),]
  exclude_components <- as.character(exclude_components$Component[!exclude_components$sig_cohort])
  
  #output
  out <- list(p, exclude_components)
  return(out)
  
}

#plot exposures grouped by parental nodes
plot_hdp_exposure_group <- function(hdpsample, group_df, incl_nonsig = T, component_names = NULL, title = NULL){
  #info:
  # group_df = data.frame dpindices,  sample names in order how it was provided to hdp and grouping
  
  # input checks
  if (!class(hdpsample) %in% c("hdpSampleChain", "hdpSampleMulti")) {
    stop("hdpsample must have class hdpSampleChain or hdpSampleMulti")
  }
  # if (!validObject(hdpsample)) stop("hdpsample not valid") # too slow on big objects
  if (length(comp_categ_counts(hdpsample)) == 0) {
    stop("No component info for hdpsample. First run hdp_extract_components")
  }
  
  dp_distn <- comp_dp_distn(hdpsample)
  ndp      <- nrow(dp_distn$mean)
  ncomp    <- ncol(dp_distn$mean)
  
  # mean exposures
  dpindices <- group_df$dpindices
  exposures <- t(dp_distn$mean[dpindices,,drop=FALSE])
  if(is.null(component_names)){
    component_names <- rownames(exposures)
  }
  rownames(exposures) <- component_names
  
  # only include significantly non-zero exposures
  if (!incl_nonsig){
    cis <- dp_distn$cred.int[dpindices]
    nonsig <- lapply(cis, function(x) which(x[1,]==0))
    for (i in 1:length(nonsig)){
      exposures[nonsig[[i]],i] <- 0
    }
    # exclude component 0
    exposures <- exposures[-1,]
  }
  
  colnames(exposures) <- group_df$sample
  
  # Number of data items per DP
  if (class(hdpsample) == "hdpSampleChain") {
    dps <- dp(final_hdpState(hdpsample))[dpindices]
    pps <- ppindex(final_hdpState(hdpsample))[dpindices]
  } else if (class(hdpsample) == "hdpSampleMulti") {
    dps <- dp(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
    pps <- ppindex(final_hdpState(chains(hdpsample)[[1]]))[dpindices]
  }
  
  numdata       <- sapply(dps, function(x) x@numdata)
  numdata_df    <- data.frame(sample = group_df$sample, numdata)
  order_samples <- numdata_df$sample[order(numdata_df$numdata, decreasing=TRUE)]
  
  
  #create data frame for plotting
  plot_data <- data.frame(reshape2::melt(exposures), type = 'Exposure')
  colnames(plot_data)[1:2] <- c('Component', 'Sample')
  plot_data <- rbind(plot_data, 
                     data.frame(Component = 'all', Sample = numdata_df$sample, value = numdata_df$numdata, type = 'Count'))
  plot_data$Sample <- factor(plot_data$Sample, levels = order_samples)
  plot_data$group  <- group_df$group[match(plot_data$Sample, group_df$sample)]
  
  component_colour <- c(brewer.pal(n = 8, 'Set3'), brewer.pal(n = 8, 'Set2'), brewer.pal(n = 8, 'Set1'))[1:ncomp]
  names(component_colour) <- component_names
  
  #plot
  ggplot(plot_data, aes(x = Sample, y = value, fill = Component)) + 
    geom_bar(stat = 'identity') +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(name = 'Component', values = c(component_colour, 'all' = '#bababa')) +
    facet_grid(type ~ group, scales = 'free', space = 'free_x') + 
    xlab('') + ylab('') +
    ggtitle(title) +
    theme_bw() + 
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
  
}



#####################################
#####           Main            #####
#####################################

#read in input_matrix
samples_input <- readRDS(input_matrix_file)
treeLayer_df  <- readRDS(treeLayer_file)

if(!is.na(nMut_cutoff)){
  samples_input <- samples_input[rowSums(samples_input) >= nMut_cutoff,]
  treeLayer_df  <- treeLayer_df[treeLayer_df$sample %in% rownames(samples_input),,drop = F]
}

#order samples_input based on treeLayer_df
samples_input <- samples_input[match(treeLayer_df$sample, rownames(samples_input)),]

#number priors
if(priors_file != 'NA'){
  prior_sigs  <- readRDS(priors_file)
  nps         <- ncol(prior_sigs) + 1
} else {
  nps <- 0
}

#dpindices
if(ncol(treeLayer_df) > 1){
  dpindex <- sum(sapply(2:ncol(treeLayer_df), function(x) length(unique(treeLayer_df[,x])))) + nps + 1
} else {
  dpindex <- nps + 1
}
dpindices <- dpindex + 1:nrow(samples_input)


# combine results from different iterations #
n_iterations <- length(grep('hdp', list.files(iteration_dir)))
chain_list <- vector("list", n_iterations)
for(i in 1:n_iterations){
  print(i)
  load(paste0(iteration_dir, 'hdp_chains_', i, '.RData'))
  chain_list[[i]] <- chains
}
mut_example_multi <- hdp_multi_chain(chain_list)


# diagnostic plots #
#--> Always remember to check that the diagnostic plots show no strong trends over the posterior samples collected!
pdf(paste0(output_dir, 'diagnosticPlots.pdf'), width = 6, height = 5)
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")
dev.off()


# Extract components (mutational signatures) #
hdpsample <- hdp_extract_components(mut_example_multi)
save(hdpsample, file = paste0(output_dir, 'hdp_results.RData'))
# load(paste0(output_dir, 'hdp_results.RData'))
# hdpsample <- hdp_results

#plot component size
pdf(paste0(output_dir, 'componentSize.pdf'), width = 7, height = 5)
par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(hdpsample, bty="L")
dev.off()

#add signatures names to priors
comp_distn <- comp_categ_distn(hdpsample)
components <- rownames(comp_distn$mean)

if(priors_file != 'NA'){
  prior_names <- grep('P', components, value = T)
  sig_index      <- as.numeric(sub('P', '', prior_names))
  sig_names      <- colnames(prior_sigs)[sig_index]
  components[grep('P', components)] <- paste0(grep('P', components, value = T), ' (', sig_names, ')')
} else {
  components[components != "0"] <- paste0('N', components[components != "0"])
}


# plot components
# SBS #
if(ncol(samples_input) == 96){
  #prepare data to plot
  sigs           <- comp_distn$mean
  colnames(sigs) <- colnames(samples_input)
  rownames(sigs) <- components
  sigs_df        <- reshape2::melt(sigs)
  colnames(sigs_df) <- c('component', 'channel', 'value')
  sigs_df$group     <- substr(as.character(sigs_df$channel), start = 3, stop = 5)
  sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels =  colnames(samples_input)),
           group = factor(group, levels = unique(sigs_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'), unique(sigs_df$group))
  strip_name_colours <- c('black','white','white','black','black','black')
  
  #plot
  pdf(paste0(output_dir, 'components.pdf'), width = 13, height = 4)
  for(x in as.character(unique(sigs_df$component))){
    plot_data <- sigs_df[sigs_df$component == x,]
    xlabels   <- paste0(substr(as.character(plot_data$channel), start = 1, stop = 1),
                        substr(as.character(plot_data$channel), start = 3, stop = 3),
                        substr(as.character(plot_data$channel), start = 7, stop = 7))
    
    p <- ggplot(plot_data, aes(x = channel, y = value, fill = group)) +
      geom_bar(stat = 'identity') + 
      facet_grid(.~ group, space = 'free_x', scales = 'free_x') +
      scale_fill_manual(name = '', values = colours, guide = 'none') +
      xlab('') +
      ylab('% SBS') +
      scale_y_continuous(expand = c(0,0,0.05,0)) +
      scale_x_discrete(labels = xlabels) +
      ggtitle(paste0('Component ', x)) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                              strip.text = element_text(face = 'bold'))
    
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
    
    grid::grid.draw(g)
    grid::grid.newpage()
  }
  dev.off()
}

# DBS #
if(ncol(samples_input) == 78){
  #prepare data to plot
  sigs           <- comp_distn$mean
  colnames(sigs) <- colnames(samples_input)
  rownames(sigs) <- components
  sigs_df        <- reshape2::melt(sigs)
  colnames(sigs_df) <- c('component', 'channel', 'value')
  sigs_df$group     <- paste0(matrix(unlist(strsplit(as.character(sigs_df$channel), '>')), ncol = 2, byrow = T)[,1], '>NN')
  sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels =  colnames(samples_input)),
           group = factor(group, levels = unique(sigs_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
                        '#e3211d', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'), unique(sigs_df$group))
  strip_name_colours <- c('black','white','black','white','black','white','black','white','black','white')
  
  #plot
  pdf(paste0(output_dir, 'components.pdf'), width = 13, height = 4)
  for(x in as.character(unique(sigs_df$component))){
    plot_data <- sigs_df[sigs_df$component == x,]
    xlabels <- matrix(unlist(strsplit(as.character(plot_data$channel), '>')), ncol = 2, byrow = T)[,2]
    
    p <- ggplot(plot_data, aes(x = channel, y = value, fill = group)) +
      geom_bar(stat = 'identity') + 
      facet_grid(.~ group, space = 'free_x', scales = 'free_x') +
      scale_fill_manual(name = '', values = colours, guide = 'none') +
      xlab('') +
      ylab('% DBS') +
      scale_y_continuous(expand = c(0,0,0.05,0)) +
      scale_x_discrete(labels = xlabels) +
      ggtitle(paste0('Component ', x)) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                              strip.text = element_text(face = 'bold'))
    
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
    
    grid::grid.draw(g)
    grid::grid.newpage()
  }
  dev.off()
}


# INDELS #
if(ncol(samples_input) == 83){
  #prepare data to plot
  sigs           <- comp_distn$mean
  colnames(sigs) <- colnames(samples_input)
  rownames(sigs) <- components
  sigs_df        <- reshape2::melt(sigs)
  colnames(sigs_df) <- c('component', 'channel', 'value')
  sigs_df$group     <- paste(matrix(unlist(strsplit(as.character(sigs_df$channel), ':')), ncol = 4, byrow = T)[,1],
                             matrix(unlist(strsplit(as.character(sigs_df$channel), ':')), ncol = 4, byrow = T)[,2],
                             matrix(unlist(strsplit(as.character(sigs_df$channel), ':')), ncol = 4, byrow = T)[,3], sep = ':')
  sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels = colnames(samples_input)),
           group = factor(group, levels = unique(sigs_df$group)),
           value = value * 100)

  #set colours
  colours <- setNames(c('#fdbe6f', '#ff8002', '#b0dd8b', '#36a12e', '#fdcab5', '#fc8a6a', '#f14432','#bc191a', '#d0e1f2', '#94c4df', '#4a98c9', '#1764ab', 
                        '#e1e1ef', '#b6b6d8', '#8683bd','#62409b'), unique(sigs_df$group))
  legend_labels <- setNames(c('1bp C Deletion', '1bp T Deletion','1bp C Insertion', '1bp T Insertion',
                              '2bp Deletion at Repeats', '3bp Deletion at Repeats', '4bp Deletion at Repeats', '5+bp Deletion at Repeats',
                              '2bp Insertion at Repeats', '3bp Insertion at Repeats', '4bp Insertion at Repeats', '5+bp Insertion at Repeats',
                              '2bp Microhomology Deletion', '3bp Microhomology Deletion', '4bp Microhomology Deletion', '5+bp Microhomology Deletion'),
                            unique(sigs_df$group))
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
  pdf(paste0(output_dir, 'components.pdf'), width = 13, height = 5)
  for(x in as.character(unique(sigs_df$component))){
    plot_data <- sigs_df[sigs_df$component == x,]
    
    p <- ggplot(plot_data, aes(x = channel, y = value, fill = group)) +
      geom_bar(stat = 'identity') + 
      facet_grid(.~group, space = 'free_x', scales = 'free_x') +
      scale_fill_manual(name = '', values = colours, labels = legend_labels) +
      xlab('Homopolymer Length / Number of Repeat Units / Microhomology Length') +
      ylab('% Indels') +
      scale_y_continuous(expand = c(0,0,0.05,0)) +
      scale_x_discrete(breaks = as.character(plot_data$channel), labels = xlabels) +
      ggtitle(paste0('Component ', x)) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15),
                              legend.position = 'bottom', strip.text = element_text(face = 'bold'))
    
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
    
    grid::grid.draw(g)
    grid::grid.newpage()
  }
  dev.off()
}



# SVs #
if(ncol(samples_input) == 32){
  #prepare data to plot
  sigs           <- comp_distn$mean
  colnames(sigs) <- colnames(samples_input)
  rownames(sigs) <- components
  sigs_df        <- reshape2::melt(sigs)
  colnames(sigs_df) <- c('component', 'channel', 'value')
  sigs_df$group <- sapply(as.character(sigs_df$channel), function(x){ unlist(strsplit(x, ':'))[1] })
  sigs_df <- sigs_df %>%
    mutate(channel = factor(channel, levels = unique(sigs_df$channel)),
           group = factor(group, levels = unique(sigs_df$group)),
           value = value * 100)
  
  #set colours
  colours <- setNames(c('#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#1f78b4', '#33a02c', '#e31a1c','#ff7f00'), unique(sigs_df$group))
  strip_name_colours <- c(rep('black', 4), rep('white', 4))
  strip_name <- matrix(unlist(strsplit(as.character(unique(sigs_df$group)), '_')), ncol = 2, byrow = T)[,2]
  
  legend_labels <- setNames(c('Clustered Deletion', 'Clustered Tandem Duplication','Clustered Inversion', 'Clustered Translocation',
                              'Non-clustered Deletion', 'Non-clustered Tandem Duplication','Non-clustered Inversion', 'Non-clustered Translocation'),
                            unique(sigs_df$group))
  
  xlabels <- as.character(sapply(unique(as.character(sigs_df$channel)), function(x){ unlist(strsplit(x, ':'))[2] }))
  xlabels <- xlabels[!is.na(xlabels)] 
  
  #plot
  pdf(paste0(output_dir, 'components.pdf'), width = 11, height = 5)
  for(x in as.character(unique(sigs_df$component))){
    plot_data <- sigs_df[sigs_df$component == x,]
    
    p <- ggplot(plot_data, aes(x = channel, y = value, fill = group)) +
      geom_bar(stat = 'identity') + 
      facet_grid(.~group, space = 'free_x', scales = 'free_x') +
      scale_fill_manual(name = '', values = colours[c(1,5,2,6,3,7,4,8)], labels = legend_labels[c(1,5,2,6,3,7,4,8)]) +
      xlab('Rearrangment size') +
      ylab('% SVs') +
      scale_y_continuous(expand = c(0,0,0.05,0)) +
      scale_x_discrete(breaks = grep('trans', as.character(plot_data$channel), invert = T, value = T), labels = xlabels) +
      ggtitle(paste0('Component ', x)) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                              strip.text = element_text(face = 'bold'), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), legend.position = 'bottom')
    
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
    grid::grid.draw(g)
    grid::grid.newpage()
  }
  dev.off()
}


# SCNA #
if(ncol(samples_input) == 48){
  #prepare data to plot
  sigs              <- comp_distn$mean
  colnames(sigs)    <- colnames(samples_input)
  rownames(sigs)    <- components
  sigs_df           <- reshape2::melt(sigs)
  colnames(sigs_df) <- c('component', 'channel', 'value')
  sigs_df$CN        <- sapply(as.character(sigs_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[1] })
  sigs_df$type      <- sapply(as.character(sigs_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[2] })
  sigs_df$length    <- sapply(as.character(sigs_df$channel), FUN = function(x){unlist(strsplit(x, ":"))[3] })
  sigs_df$value     <- sigs_df$value * 100
  sigs_df$group     <- paste0(sigs_df$CN, ":", sigs_df$type)
  sigs_df$group     <- factor(sigs_df$group, levels = unique(sigs_df$group))
  
  sigs_df$type      <- factor(sigs_df$type, levels = c("homdel", "LOH", "het"))
  sigs_df$length    <- as.character(sigs_df$length)
  sigs_df$channel   <- factor(sigs_df$channel, levels = unique(sigs_df$channel))
  
  #set colours
  colours       <- setNames(c('#fb8072', '#ffd92f','#66c2a5', '#e78ac3',  '#8da0cb', '#fc8d62', '#1b9e77', '#e7298a',  '#7570b3', '#d95f02'), unique(sigs_df$group))
  legend_labels <- setNames(c('homozygous deletion', 'LOH & CN1', 'LOH & CN2', 'LOH & CN3-4', 'LOH & CN5-8', 'LOH & CN9+', 
                              'heterozygous & CN2', 'heterozygous & CN3-4', 'heterozygous & CN5-8', 'heterozygous & CN9+'), unique(sigs_df$group))
  
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
  pdf(paste0(output_dir, 'components.pdf'), width = 11, height = 5)
  for(x in as.character(unique(sigs_df$component))){ 
    plot_data <- sigs_df[sigs_df$component == x,]
    
    p <- ggplot(plot_data, aes(x = channel, y = value, fill = group)) +
      geom_bar(stat = 'identity') + 
      facet_grid(.~group, space = 'free_x', scales = 'free_x') +
      scale_fill_manual(name = '', values = colours[c(1,2,3,7,4,8,5,9,6,10)], labels = legend_labels[c(1,2,3,7,4,8,5,9,6,10)]) +
      xlab('Segment Length') +
      ylab('% CN segments') +
      scale_y_continuous(expand = c(0,0,0.05,0)) +
      scale_x_discrete(breaks = signatures_df$channel, labels = xlabels) +
      ggtitle(paste0('Component ', x)) +
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
    grid::grid.draw(g)
    grid::grid.newpage()
  }
  dev.off()
}


# plot exposures
pdf(paste0(output_dir, 'componentExposures.pdf'), width = 15, height = 7)
plot_dp_comp_exposure(hdpsample, main_text="Exposures (excluding non-significant)",
                      dpindices=dpindices,
                      col= c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(12, "Set3")),
                      incl_nonsig = F,
                      incl_numdata_plot = T,
                      ylab_numdata = 'Count', ylab_exp = 'Exposure',
                      leg.title = 'Components')
plot_dp_comp_exposure(hdpsample, main_text="Exposures (including non-significant)",
                      dpindices=dpindices,
                      col= c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(12, "Set3")),
                      incl_nonsig = T,
                      incl_numdata_plot = T,
                      ylab_numdata = 'Count', ylab_exp = 'Exposure',
                      leg.title = 'Components')
dev.off()



# seperating samples with different parental nodes (e.g. subtypes) if used
if(ncol(treeLayer_df) > 1){
  p_excl_nonsig <- lapply(2:ncol(treeLayer_df), function(x){
    if(length(unique(treeLayer_df[,x])) > 20){ return(NULL) }
    group_df <- data.frame(dpindices = dpindices, sample = treeLayer_df$sample, group = treeLayer_df[,x])
    plot_hdp_exposure_group(hdpsample, group_df = group_df, incl_nonsig = F, 
                            component_names = components, title = "Exposures (excluding non-significant)")
  })
  
  p_incl_nonsig <- lapply(2:ncol(treeLayer_df), function(x){
    if(length(unique(treeLayer_df[,x])) > 20){ return(NULL) }
    group_df <- data.frame(dpindices = dpindices, sample = treeLayer_df$sample, group = treeLayer_df[,x])
    plot_hdp_exposure_group(hdpsample, group_df = group_df, incl_nonsig = T, 
                            component_names = components, title = "Exposures (including non-significant)")
  })
  
  pdf(paste0(output_dir, 'componentExposures_perGroup.pdf'), width = 17, height = 7)
  print(p_excl_nonsig)
  print(p_incl_nonsig)
  dev.off()
}


# plot exposures as boxplot
qc_activity <- plot_hdp_exposure_boxplot(hdpsample, dpindices, component_names = components, sig_active_cutoff = sigActivity_cutoff, cohort_threshold = cohort_cutoff)
exclude_components <- qc_activity[[2]]
exclude_components <- as.character(sapply(exclude_components, function(x){unlist(strsplit(x, ' '))[1]}))

pdf(paste0(output_dir, 'componentExposures_boxplot.pdf'), width = 10, height = 4)
qc_activity[[1]]
dev.off()


# reconstruction error
comp_distn <- comp_categ_distn(hdpsample)
signatures <- comp_distn$mean
dp_distn   <- comp_dp_distn(hdpsample)
exposures  <- dp_distn$mean[dpindices,]

observed_mutLoad <- samples_input
expected_mutLoad <- exposures %*% signatures
rownames(expected_mutLoad) <- rownames(observed_mutLoad)

RMSE  <- sapply(rownames(observed_mutLoad), function(x) rmse(as.numeric(observed_mutLoad[x,]), as.numeric(expected_mutLoad[x,])))
nRMSE <- RMSE / rowMeans(observed_mutLoad)
cosineSimilarity <- sapply(rownames(observed_mutLoad), function(x) cosine(as.numeric(observed_mutLoad[x,]), as.numeric(expected_mutLoad[x,])))

save(cosineSimilarity, RMSE, nRMSE, file = paste0(output_dir, "reconstructionError.RData"))

#plot reconstruction error
plot_data <- rbind(data.frame(sample = names(cosineSimilarity), value = as.numeric(cosineSimilarity), type = 'cosineSimilarity'),
                   data.frame(sample = names(nRMSE), value = as.numeric(nRMSE), type = 'nRMSE'),
                   data.frame(sample = names(RMSE), value = as.numeric(RMSE), type = 'RMSE'))

order_samples <- data.frame(sample = names(cosineSimilarity), value = as.numeric(cosineSimilarity), type = 'cosineSimilarity')
order_samples <- order_samples[order(order_samples$value, decreasing = T),]
plot_data$sample <- factor(plot_data$sample, levels = order_samples$sample)

p <- ggplot(plot_data, aes(x = sample, y = value, fill = type)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#fdae61', '#f46d43', '#a50026')) +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(type ~ ., scales = 'free_y') +
  xlab('Patients') + ylab('') +
  ggtitle('Reconstruction Error per Patient') +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     legend.position = 'none')
pdf(paste0(output_dir, "reconstructionError.pdf"), width = 10, height = 5)
plot(p)
dev.off()





####### save results #######
hdp_results <- hdpsample
save(hdp_results, file = paste0(output_dir, 'hdp_results.RData'))


# save signatures / components ##
comp_distn <- comp_categ_distn(hdpsample)
signatures <- comp_distn$mean
colnames(signatures) <- colnames(samples_input)
signatures <- t(signatures)

write.table(signatures, file = paste0(output_dir, 'components.txt'), sep ="\t")

#exclude Component 0 and components with low activity
signatures <- signatures[,!colnames(signatures) %in% c('0', exclude_components)]

if(priors_file != 'NA'){
  prior_names <- grep('P', colnames(signatures), value = T)
  sig_index      <- as.numeric(sub('P', '', prior_names))
  SBS_names      <- colnames(prior_sigs)[sig_index]
  colnames(signatures)[grep('P', colnames(signatures))] <- SBS_names
} else {
  colnames(signatures) <- paste0('N', colnames(signatures))
}

write.table(signatures, file = paste0(output_dir, 'signatures.txt'), sep ="\t")

## save exposures ##
dp_distn  <- comp_dp_distn(hdpsample)
exposures <- dp_distn$mean[dpindices,]
rownames(exposures) <- rownames(samples_input)

write.table(exposures, file = paste0(output_dir, 'componentExposure.txt'), sep ="\t")

mutBurden <- data.frame(sample = names(rowSums(samples_input)), nMut = as.numeric(rowSums(samples_input)))
mutBurden <- mutBurden[match(rownames(exposures), mutBurden$sample),]
exposures_counts <- round(exposures * mutBurden$nMut)

write.table(exposures_counts, file = paste0(output_dir, 'componentExposure_counts.txt'), sep ="\t")


#exclude Component 0 and components with low activity
# exclude non-significant exposures
cis <- dp_distn$cred.int[dpindices]
nonsig <- lapply(cis, function(x) which(x[1,]==0))
signf_exposures <- exposures
for (i in 1:length(nonsig)){
  signf_exposures[i, nonsig[[i]]] <- 0
}

exposures       <- exposures[,!colnames(exposures) %in% c('0', exclude_components)]
signf_exposures <- signf_exposures[,!colnames(signf_exposures) %in% c('0', exclude_components)]

if(priors_file != 'NA'){
  prior_names <- grep('P', colnames(exposures), value = T)
  sig_index      <- as.numeric(sub('P', '', prior_names))
  SBS_names      <- colnames(prior_sigs)[sig_index]
  colnames(exposures)[grep('P', colnames(exposures))] <- SBS_names
  colnames(signf_exposures)[grep('P', colnames(signf_exposures))] <- SBS_names
} else {
  colnames(exposures) <- paste0('N', colnames(exposures))
  colnames(signf_exposures) <- paste0('N', colnames(signf_exposures))
}

write.table(exposures, file = paste0(output_dir, 'signatureExposures.txt'), sep ="\t")
write.table(signf_exposures, file = paste0(output_dir, 'signatureExposures_signif.txt'), sep ="\t")

mutBurden <- data.frame(sample = names(rowSums(samples_input)), nMut = as.numeric(rowSums(samples_input)))
mutBurden <- mutBurden[match(rownames(exposures), mutBurden$sample),]
exposures_counts       <- round(exposures * mutBurden$nMut)
signf_exposures_counts <- round(signf_exposures * mutBurden$nMut)

write.table(exposures_counts, file = paste0(output_dir, 'signatureExposures_counts.txt'), sep ="\t")
write.table(signf_exposures_counts, file = paste0(output_dir, 'signatureExposures_signif_counts.txt'), sep ="\t")
  









