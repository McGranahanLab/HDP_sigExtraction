#!/bin/bash 

###########################################################
######        HDP Signature Extraction Wrapper       ######
###########################################################

### paramters ###
outputBase="."                                     #create set output directory which need to include a folder called input_files that contains the mutation count matrices and treeLayer files from the setup file
setUpFile='runSetup.txt'                           #file located in outputBase/input_files/ that contains set up for which input matrix should be run with which treeLayer file
refFile="./data/COSMIC_v3.2_SBS_GRCh37.txt"        #reference file to compare the extracted signatures to --> needs to be changed when other types of signatures than SBS are analysed
normFile="./data/trinucContext_exome_hg19.rds"     #if normalisation of  signatures should be applied the a normFile needs to be provided, otherwise set NA
priorFile="NA"                                     #if priors included in the run the a priorFile needs to be provided, otherwise set NA
nMutCutoff=50                                      #samples with less mutations that specified in this parameter will be excluded from the analysis

# signature activity
sigActivityThreshold=0                             #minimum signature exposure cutoff
cohortThreshold=0                                  #only signatures that are active in a higher proportion of samples in the cohort than specified in this threshold will be used for further analyses

# HDP sampling
burnin=10000
nPosterior=100
space=200
cpiter=3
iterations=15

# signature comparison
cosineSimThreshold=0.9
maxiterEM=1000
EMfracThreshold=0.1

#script directory
scriptDir='./R/'


#####################
####     Main    ####
#####################

while read l;
	do
		lineArray=( $l )
		countMatrix=${lineArray[0]} 
		treeLayers=${lineArray[1]} 

		#input files
		input_96matrix_file="$outputBase/input_files/$countMatrix"
		treeLayer_file="$outputBase/input_files/$treeLayers.rds"

		#cancer signature file
		cancerSigsFile="$outputBase/input_files/lungSigs_list.rds"


		### create outputDir based on input parameters ###
		priorsName="with_priors"
		if [ ${priorFile} == "NA" ]
		then 
			priorsName="without_priors"
		fi

		nMutName='without_minMut'
		if [ ${nMutCutoff} != "NA" ]
		then
			nMutName="minMut_$nMutCutoff"
		fi

		hdpParamsName="burnin${burnin}_n${nPosterior}_space${space}_cpiter${cpiter}"

		outputDir="$outputBase/$treeLayers/$priorsName/$nMutName/$hdpParamsName/"
		iterationDir="$outputDir/iterations/"
		logDir="$outputDir/logs/"

		if [ ! -d "$iterationDir" ]
		then 
			mkdir -p $iterationDir
		fi

		if [ ! -d "$logDir" ]
		then 
			mkdir $logDir
		fi


		### run hdp sampling ###
		stdout="$logDir/hdp_chainSampling.stdout"
		stderr="$logDir/hdp_chainSampling.stderr"

		jobid1=$(sbatch --parsable --time=3-00:00:00 -c 2 --array=1-$iterations --mem 12G -J HDPchain -o $stdout -e $stderr $scriptDir/run_HDP_chainSampling.R ${input_96matrix_file} $iterationDir ${treeLayer_file} $priorFile $nMutCutoff $burnin $nPosterior $space $cpiter)


		### process hdp results ###
		stdout="$logDir/hdp_processing.stdout"
		stderr="$logDir/hdp_processing.stderr"

		jobid2=$(sbatch --parsable --dependency="afterok:$jobid1" --time=1-00:00:00 -c 4 --mem 12G -J HDPprocess -o $stdout -e $stderr $scriptDir/run_HDP_processing.R ${input_96matrix_file} $outputDir ${treeLayer_file} $priorFile $nMutCutoff $sigActivityThreshold $cohortThreshold)

		### normalise signatures if required ###
		if [ ${normFile} != "NA" ]
		then
			stdout="$logDir/hdp_sigNormalising.stdout"
		  stderr="$logDir/hdp_sigNormalising.stderr"

		  jobid3=$(sbatch --parsable --dependency="afterok:$jobid2" --time=01:00:00 -c 2 --mem 8G -J HDPnorm -o $stdout -e $stderr $scriptDir/run_HDP_sigNormalising.R $outputDir $normFile) 
		    

      stdout="$logDir/hdp_normSigs_comparingCOSMIC.stdout"
			stderr="$logDir/hdp_normSigs_comparingCOSMIC.stderr"
			  
			jobid4=$(sbatch --parsable --dependency="afterok:$jobid3" --time=1-00:00:00 -c 2 --mem 8G -J HDPcompare -o $stdout -e $stderr $scriptDir/run_HDP_comparing.R $outputDir $cosineSimThreshold $maxiterEM $EMfracThreshold $refFile "NA" "TRUE")
		fi



		### compare hdp to comsic ###
		stdout="$logDir/hdp_comparingCOSMIC.stdout"
		stderr="$logDir/hdp_comparingCOSMIC.stderr"
			  
		jobid5=$(sbatch --parsable --dependency="afterok:$jobid2" --time=1-00:00:00 -c 2 --mem 8G -J HDPcompare -o $stdout -e $stderr $scriptDir/run_HDP_comparing.R $outputDir $cosineSimThreshold $maxiterEM $EMfracThreshold $refFile "NA" "FALSE")
		

		stdout="$logDir/hdp_comparingCOSMIC_cancerSigs.stdout"
		stderr="$logDir/hdp_comparingCOSMIC_cancerSigs.stderr"

		jobid6=$(sbatch --parsable --dependency="afterok:$jobid5" --time=1-00:00:00 -c 2 --mem 8G -J HDPcompare -o $stdout -e $stderr $scriptDir/run_HDP_comparing.R $outputDir $cosineSimThreshold $maxiterEM $EMfracThreshold $refFile $cancerSigsFile "FALSE")
		

	done < $outputBase/input_files/$setUpFile



