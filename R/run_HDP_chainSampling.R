#!/usr/bin/env Rscript

#############################################################################################
######        Run multiple posterior sampling chains for HDP  signature analysis       ######
#############################################################################################
#--> submit script  multiple times with wrapper to run in parallel 

#options and libraries
options(stringsAsFactors = F)
library(hdp)


#parameters
cmdArgs <- commandArgs(trailingOnly = TRUE)
input_matrix_file     <- cmdArgs[1]
output_dir            <- cmdArgs[2]
treeLayer_file        <- cmdArgs[3]
priors_file           <- cmdArgs[4]             #NA = no priors
norm_file             <- cmdArgs[5]             #NA = no normalisation needed 
nMut_cutoff           <- as.numeric(cmdArgs[5]) #NA = no minimum mutations cutoff
burnin                <- as.numeric(cmdArgs[6])
n_posterior           <- as.numeric(cmdArgs[7])
space                 <- as.numeric(cmdArgs[8])
cpiter                <- as.numeric(cmdArgs[9])
iter                  <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))




#####################################
#####           Main            #####
#####################################

#read in input and treeLayers
samples_input <- readRDS(input_matrix_file)
treeLayer_df  <- readRDS(treeLayer_file)

if(!is.na(nMut_cutoff)){
  samples_input <- samples_input[rowSums(samples_input) >= nMut_cutoff,]
  treeLayer_df  <- treeLayer_df[treeLayer_df$sample %in% rownames(samples_input),,drop = F]
}

#order samples_input based on treeLayer_df
samples_input <- samples_input[match(treeLayer_df$sample, rownames(samples_input)),]



##############
###  initialise HDP structure  
#############

### with priors ###
if(priors_file != 'NA'){
  prior_sigs  <- readRDS(priors_file)
  nps         <- ncol(prior_sigs)
  
  #signature priors
  ttype_prior <- hdp_prior_init(prior_distn = prior_sigs,
                                prior_pseudoc = rep(1000, nps), 
                                hh=rep(1, ncol(samples_input)),
                                alphaa=c(1, 1),
                                alphab=c(1, 1))
  
  if(ncol(treeLayer_df) > 1){
    nconParam   <- sum(sapply(2:ncol(treeLayer_df), function(x) length(unique(treeLayer_df[,x])))) + nps + 3
  } else {
    nconParam   <- nps + 3
  }
  
  ttype_prior <- hdp_addconparam(ttype_prior,
                                 alphaa = rep(1, nconParam),
                                 alphab = rep(1, nconParam))
  
  #add nodes for layes in treeLayer_df
  for(i in rev(1:ncol(treeLayer_df))){
    if(i == ncol(treeLayer_df)){
      numdp  <- 1 + length(unique(treeLayer_df[,i]))
      ppindex <- c(1, rep(1+nps+1, length(unique(treeLayer_df[,i]))))
      cpindex <- c(3, rep(4, length(unique(treeLayer_df[,i]))))
    } else {
      numdp    <- length(unique(treeLayer_df[,i]))
      index_df <- treeLayer_df[!duplicated(treeLayer_df[,i]),,drop = F]
      index    <- as.numeric(factor(index_df[,i+1], levels = unique(index_df[,i+1])))
      ppindex  <- index + max(ttype_prior@ppindex)
      cpindex  <- index + max(ttype_prior@cpindex)
    }
    
    ttype_prior <- hdp_adddp(ttype_prior,
                             numdp = numdp,                                                                                                                                 
                             ppindex = ppindex,
                             cpindex = cpindex)
    
    }
 
  
  #add data
  dpindex     <- (length(ttype_prior@ppindex) - nrow(samples_input))
  ttype_prior <- hdp_setdata(ttype_prior,
                             dpindex = dpindex + 1:nrow(samples_input),
                             samples_input)
  
  
  #active added nodes
  nLayers <- sum(sapply(1:ncol(treeLayer_df), function(x) length(unique(treeLayer_df[,x]))))
  ttype_activated <- dp_activate(ttype_prior,
                                 dpindex = (1+nps+1)+0:nLayers,
                                 initcc = nps+10,
                                 seed = iter*1000)
  
  
}



### without priors ###
if(priors_file == 'NA'){
  
  ppindex <- 0
  cpindex <- 1
  for(i in rev(1:ncol(treeLayer_df))){
    if(i == ncol(treeLayer_df)){
      ppindex <- c(ppindex, rep(1, length(unique(treeLayer_df[,i]))))
      cpindex <- c(cpindex, rep(3, length(unique(treeLayer_df[,i]))))
    } else{
      index_df <- treeLayer_df[!duplicated(treeLayer_df[,i]),,drop = F]
      index    <- as.numeric(factor(index_df[,i+1], levels = unique(index_df[,i+1])))
      ppindex  <- c(ppindex, index + max(ppindex))
      cpindex  <- c(cpindex, index + max(cpindex))
    }
  }

  # initialise 
  hdp_mut <- hdp_init(ppindex = ppindex,
                      cpindex = cpindex,
                      hh=rep(1, ncol(samples_input)), 
                      alphaa=rep(1, max(cpindex)), 
                      alphab=rep(1, max(cpindex))) 
  
  # add data to leaf nodes (one per sample, in row order of mut_count)
  dpindex     <- (length(hdp_mut@ppindex) - nrow(samples_input)) + 1
  hdp_mut <- hdp_setdata(hdp_mut, 
                         dpindex = dpindex:numdp(hdp_mut), 
                         samples_input) 
  
  #active added nodes
  ttype_activated <- dp_activate(hdp_mut,
                                 dpindex = 1:numdp(hdp_mut),
                                 initcc = 10,
                                 seed = iter*1000)
  
}



############
### run posterior sampling chains
############

#run chain sampling
chains <- hdp_posterior(ttype_activated,
                        burnin = burnin,
                        n = n_posterior,
                        space = space,
                        cpiter = cpiter,
                        seed = iter*1e6,
                        verbosity = 1)


#save results and combine multiple runs in different script
save(chains, file = paste0(output_dir, '/hdp_chains_', iter, '.RData'))




