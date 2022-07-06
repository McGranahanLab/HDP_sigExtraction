#!/usr/bin/env Rscript

###########################################################################################################
######                Normalise signatures based on trinucleotide context background                 ######
###########################################################################################################

# options and libraries
options(stringsAsFactors = F)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)



#parameters
cmdArgs               <- commandArgs(trailingOnly = TRUE)
output_dir            <- cmdArgs[1]          
normalisation_file    <- cmdArgs[2]            


###################################
#####       Functions         #####
###################################

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




#####################################
#####           Main            #####
#####################################

#create norm folder
norm_dir <- paste0(output_dir, '/normSignatures/')
if(!file.exists(norm_dir)){
  dir.create(norm_dir)
}

#read in signatures and components
components <- read.table(paste0(output_dir, 'components.txt'), header = T, sep = '\t')
colnames(components) <- sub('X', '', colnames(components))
signatures <- read.table(paste0(output_dir, 'signatures.txt'), header = T, sep = '\t')

#load trinuc context count
trinucCounts <- readRDS(normalisation_file)
trinucCounts$ratio <- trinucCounts$genome_count / trinucCounts$region_count



#assign trinucContext to trinucContextMutations
trinuc_df <- data.frame(trinuc_mut = rownames(signatures),
                        base1 = substr(rownames(signatures), start = 1, stop = 1),
                        base2 = substr(rownames(signatures), start = 3, stop = 3),
                        base3 = substr(rownames(signatures), start = 7, stop = 7))
trinuc_df$trinuc_context <- paste0(trinuc_df$base1, trinuc_df$base2, trinuc_df$base3)
trinuc_df <- trinuc_df %>% left_join(trinucCounts, by = 'trinuc_context')



#######   normalise signatures   #######
norm_signatures <- signatures * trinuc_df$ratio[match(rownames(signatures), trinuc_df$trinuc_mut)]
norm_signatures <- t(norm_signatures) / colSums(norm_signatures)
norm_signatures <- t(norm_signatures)
write.table(norm_signatures, file = paste0(norm_dir, 'normSignatures.txt'), sep ="\t")

#plot
plot_data <- reshape2::melt(data.frame(channel = rownames(norm_signatures), norm_signatures))
pdf(paste0(norm_dir, 'normSignatures.pdf'), width = 12, height = 4)
for(x in unique(plot_data$variable)){
  sub <- plot_data[plot_data$variable == x, c('channel', 'value')]
  grid::grid.draw(plot_SBS(sub, title = paste0('Signature ', x)))
  grid::grid.newpage()
}
dev.off()




#######   normalise components    #######
norm_components <- components * trinuc_df$ratio[match(rownames(components), trinuc_df$trinuc_mut)]
norm_components <- t(norm_components) / colSums(norm_components)
norm_components <- t(norm_components)
write.table(norm_components, file = paste0(norm_dir, 'normComponents.txt'), sep ="\t")


#plot
plot_data <- reshape2::melt(data.frame(channel = rownames(norm_components), norm_components))
plot_data$variable <- sub('X0', '0', plot_data$variable)
pdf(paste0(norm_dir, 'normComponents.pdf'), width = 12, height = 4)
for(x in unique(plot_data$variable)){
  sub <- plot_data[plot_data$variable == x, c('channel', 'value')]
  grid::grid.draw(plot_SBS(sub, title = paste0('Component ', x)))
  grid::grid.newpage()
}
dev.off()


  

  
  
  
  
