# HDP_sigExtraction
Pipeline to extract de novo mutational signatures using HDP.


This pipeline is based on the R-package hdp developed by Nicola Roberts (https://github.com/nicolaroberts/hdp). Documentation about how to apply hdp for signature extraction is provided by  the correpsonding vingette (https://github.com/nicolaroberts/hdp/blob/master/vignettes/mutation_signatures.Rmd).

## Getting started & Documentation
1) Install the R-packge hdp from github like described on the website https://github.com/nicolaroberts/hdp

2) The following three input files are required to run the pipeline:
* input matrix with the mutation counts (e.g. 96-trinucleotide counts for SBS) for each sample, with the mutation categories as columns and the samples as rows. An example can be found in examples input_96matrix_tx100.rds which contains the 96-trinucleotide counts for the TARCERx 100 cohort.
* treeLayer file which is a data frame that contains the sample names in the first column and any number of categorical variables should be considered in the dependency tree. One example with histology as tree dependecy can be found in examples treeLayers_cancerType_tx100.rds
* setUp file that includes the information about what input matrix should be run together with what treeLayer file. This is very helpful in cases where multiple different treeLayers should be applied for the same cohort of patients.

3) Optional input files:
* A file with a subset of signatures that should be used as priors.
* A file with background information, e.g. trinucleotide counts across the genome, that should be used for normalisation of the signatures
* A file with reference signatures that will be used to compare extracted signatures to. The COSMIC v3.2 signature file is set as default in the wrapper script.

All input files need to be saved in a folder called input_files in your output directory. For the results the pipeline will create a subfolder with a name that states the different parameters that were used to run this pipeline.

4) Adaptation of the HDP_wrapper.sh script:
* If files were moved around, make sure all paths were updated.
* This script contains different paremeters specific for the signature extraction and comparisons which can be updated by the user.
        
        Signature activity thresholds:
        sigActivityThreshold = 0.05
        cohortThreshold = 0.01

        HDP sampling parameters:
        burnin=10000
        nPosterior=100
        space=200
        cpiter=3
        iterations=15

        Signature comparison parameters:
        cosineSimThreshold=0.9
        maxiterEM=1000
        EMfracThreshold=0.1
        
5) Run the HDP_wrapper.sh script to run the pipeline.

