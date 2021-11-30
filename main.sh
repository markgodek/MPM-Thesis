#!/bin/bash

## ran using
# The Genome Analysis Toolkit (GATK) v4.1.8.0
# HTSJDK Version: 2.22.0
# Picard Version: 2.22.8

#DEBUG="True"
DEBUG="False"

#RESELECT="Somatic"
#RESELECT="Germline"
#RESELECT="Both"
RESELECT="False"

ONESCRIPT="True"
#ONESCRIPT="False"

## load modules

module load fastqc trimmomatic/0.36 bwa/0.7.16a samtools/1.4.1 bcftools bedtools vcftools

## static directories

SCRIPT_DIR=user_dir/example_familial/scripts/small_pipeline
EXOME_DIR=/mntraw_data_path/exome                         # a directory which contains all raw fastq reads to be input with filename pattern "HITS######_1_1.gne.fastq.gz"
#OUT=lab_path/oct26_run                   # the root directory which the completed files will be output to
OUT=lab_path/variant_calling_pipeline
AIMS_DIR=lab_path/AimIntervals                # file which contains GATK or picard-style intervals files of genes of interest

## Files which must be provided

SAM_PATIENT_LINK=raw_data_path/SAM_patient_link.csv       # CSV file with one entry per line. File should have a header and columns are patient ID, tumor sample ID, normal sample ID
REF=/mntref_genome_dir/human/BWA-MEM/b37_decoys/hs37d5.fa
INTERVALS=/mntref_genome_dir/human/hs37d5_interval_files/hs37d5.genome.exons.gatk-style.intervals
NO_DECOY_INTERVALS=/mntref_genome_dir/human/hs37d5_interval_files/hs37d5.genome.exons.no-decoys.gatk-style.intervals
RENAME_CHRs=lab_path/Annotation_data_sources/hg19tob37.chr_name.map

#known sites - One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.
#Note: the pipeline will need to be modified if number of known sites will be increased or decreased as Hg19 ref genome falls out of use for SNP calling
#Note: to generate GERMLINE_RESOURCE_BIALLELIC_COMMON_SNPS run this command on gnomad allele frequency raw sites
#gatk SelectVariants  --output af-only-gnomad.biallelic.vcf.gz --restrict-alleles-to BIALLELIC --variant /mntref_genome_dir/human/known_sites/af-only-gnomad.raw.sites.b37.vcf.gz

COMMON_SNPS=ref_genome_dir/human/known_sites/GRCh37.p13/common_all_20180423.vcf.gz
ONEK_SNPS=ref_genome_dir/human/known_sites/broad_bundle/1000G_phase1.snps.high_confidence.b37.vcf.gz 	# likely consists of all SNPs from 1000 Genomes project, filters too much
DBSNP=/mntref_genome_dir/human/known_sites/broad_bundle/dbsnp_138.b37.vcf.gz
ONEK_INDELS=/mntref_genome_dir/human/known_sites/broad_bundle/1000G_phase1.indels.b37.vcf
MILLS_INDELS=/mntref_genome_dir/human/known_sites/broad_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf
GERMLINE_RESOURCE_BIALLELIC_COMMON_SNPS=/mntref_genome_dir/human/known_sites/af-only-gnomad.biallelic.vcf.gz
OMNI=/mntref_genome_dir/human/known_sites/broad_bundle/1000G_omni2.5.b37.vcf
HAPMAP=/mntref_genome_dir/human/known_sites/broad_bundle/hapmap_3.3.b37.vcf
ANNOTATIONS=/mntlab_path/Annotation_data_sources/funcotator_dataSources.v1.6.20190124s

GERMLINE_RESOURCE=/mntref_genome_dir/human/known_sites/af-only-gnomad.raw.sites.b37.vcf.gz
NORM_GERMLINE_RESOURCE=/mntref_genome_dir/human/known_sites/norm.af-only-gnomad.raw.sites.b37.vcf.gz

if [[ "${DEBUG}" == "True" ]]
then 
	echo "STARTING PIPELINE IN DEBUG MODE"

	OUT=lab_path/debug_small_pipeline
	EXOME_DIR=/mntraw_data_path/exome/small                   # directory with only first million FASTQ reads, only for patients in debug_link.csv to speed up debugging
	#EXOME_DIR=/mntraw_data_path/exome/supersmall
	SAM_PATIENT_LINK=raw_data_path/debug_link.csv            # CSV link with smaller dataset for debugging
	
	#rm -r ${OUT}
fi

sh ${SCRIPT_DIR}/patient_parser.sh ${SAM_PATIENT_LINK} ${SCRIPT_DIR}            # parse the link file to get create files that let us count samples and run scripts based on indicies

## Counts derived from files

SAMPLE_COUNT=$(wc -l ${SCRIPT_DIR}/ALL_SAMPLES.tmp | awk '{print $1}')
PATIENT_COUNT=$(wc -l ${SCRIPT_DIR}/PATIENT_LIST.tmp | awk '{print $1}')
TUMOR_COUNT=$(wc -l ${SCRIPT_DIR}/TUMOR_ONLY_LINK.csv | awk '{print $1-1}')
NORMAL_COUNT=$(wc -l ${SCRIPT_DIR}/NORMAL_SAMPLES.tmp | awk '{print $1}')

LOGS=${OUT}/job_results
PROC_DIR=${OUT}/processing
som_ponDB=${OUT}/som_ponDB
germ_ponDB=${OUT}/germ_ponDB
PON_VCFs=${OUT}/pon
SOMATIC=${OUT}/somatic
GERMLINE=${OUT}/germline
AIMS_OUT=${OUT}/Aims
AF_DIR=${OUT}/AF

mkdir -p ${LOGS} ${PROC_DIR} ${PON_VCFs} ${SOMATIC} ${GERMLINE} ${AIMS_OUT} ${AF_DIR}

if [[ "${ONESCRIPT}" == "False" ]] && [[ "${RESELECT}" == "False" ]]
then

jid_process=$(bsub \
-n 6 \
-q big \
-W 336:00 \
-o ${LOGS}/process_%J.out \
-e ${LOGS}/process_%J.err \
-J "process[1-${SAMPLE_COUNT}]" ${SCRIPT_DIR}/process.lsf \
${EXOME_DIR} ${REF} ${INTERVALS} ${SCRIPT_DIR}/ALL_SAMPLES.tmp ${PROC_DIR} ${ONEK_INDELS} \
${MILLS_INDELS} ${DBSNP} ${SCRIPT_DIR}/NORMAL_SAMPLES.tmp ${PON_VCFs} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_somDB=$(bsub \
-n 4 \
-q big \
-W 336:00 \
-o ${LOGS}/somDB_%J.out \
-e ${LOGS}/somDB_%J.err \
-w "done(${jid_process[*]})" \
-J "somDB[1]" ${SCRIPT_DIR}/genomicsDBimport.lsf \
${PON_VCFs} ${SCRIPT_DIR}/NORMAL_SAMPLES.tmp ${som_ponDB} ${REF} \
${INTERVALS} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_somatic=$(bsub \
-n 1 \
-q big \
-W 334:00 \
-o ${LOGS}/som_analysis_%J.out \
-e ${LOGS}/som_analysis_%J.err \
-w "done(${jid_somDB[*]})" \
-J "som_var[1-${TUMOR_COUNT}]" ${SCRIPT_DIR}/analysis.lsf \
${SCRIPT_DIR}/TUMOR_ONLY_LINK.csv ${PROC_DIR} ${SOMATIC} ${REF} ${INTERVALS} \
${som_ponDB} ${GERMLINE_RESOURCE} ${PON_VCFs} ${GERMLINE_RESOURCE_BIALLELIC_COMMON_SNPS} \
${ANNOTATIONS} ${RENAME_CHRs} ${AIMS_DIR} ${AIMS_OUT} \
${DBSNP} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_haplo=$(bsub \
-n 1 \
-q big \
-W 336:00 \
-o ${LOGS}/haplo_%J.out \
-e ${LOGS}/haplo_%J.err \
-w "done(${jid_process[*]})" \
-J "haplo[1-${NORMAL_COUNT}]" ${SCRIPT_DIR}/haplotype.lsf \
${PROC_DIR} ${SCRIPT_DIR}/NORMAL_SAMPLES.tmp ${GERMLINE} \
${REF} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_germDB=$(bsub \
-n 6 \
-q big-multi  \
-W 336:00 \
-o ${LOGS}/germdb_%J.out \
-e ${LOGS}/germdb_%J.err \
-w "done(${jid_haplo[*]})" \
-J "germDB[1]" ${SCRIPT_DIR}/genomicsDBimport.lsf \
${GERMLINE} ${SCRIPT_DIR}/NORMAL_SAMPLES.tmp ${germ_ponDB} ${REF} \
${INTERVALS} True | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_jointcall=$(bsub \
-n 4 \
-q big-multi \
-W 336:00 \
-o ${LOGS}/jointcall_%J.out \
-e ${LOGS}/jointcall_%J.err \
-w "done(${jid_germDB[*]})" \
-J "jointcall[1]" ${SCRIPT_DIR}/jointcall-cohort.lsf \
${REF} ${INTERVALS} ${GERMLINE} ${germ_ponDB} ${HAPMAP} ${OMNI} ${ONEK_SNPS} \
${DBSNP} ${ANNOTATIONS} ${RENAME_CHRs} ${AIMS_DIR} \
${AIMS_OUT} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_reselect_som=$(bsub \
-n 4 \
-q big-multi \
-W 24:00 \
-o ${LOGS}/reselect_som_%J.out \
-e ${LOGS}/reselect_som_%J.err \
-w "done(${jid_somatic[*]})" \
-J "reselect_som[1-${TUMOR_COUNT}]" ${SCRIPT_DIR}/reselect.lsf \
${REF} ${AIMS_DIR} ${AIMS_OUT} ${SCRIPT_DIR}/TUMOR_PATIENTS.tmp \
${COMMON_SNPS} ${DBSNP} ${SOMATIC} Somatic ${AF_DIR} \
${NORM_GERMLINE_RESOURCE} ${SCRIPT_DIR} ${NO_DECOY_INTERVALS} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_reselect_germ=$(bsub \
-n 4 \
-q big-multi \
-W 24:00 \
-o ${LOGS}/reselect_germ_%J.out \
-e ${LOGS}/reselect_germ_%J.err \
-w "done(${jid_jointcall[*]})" \
-J "reselect_germ[1]" ${SCRIPT_DIR}/reselect.lsf \
${REF} ${AIMS_DIR} ${AIMS_OUT} ${SCRIPT_DIR}/TUMOR_PATIENTS.tmp \
${COMMON_SNPS} ${DBSNP} ${GERMLINE} Germline ${AF_DIR} \
${NORM_GERMLINE_RESOURCE} ${SCRIPT_DIR} ${NO_DECOY_INTERVALS} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

jid_cleanup=$(bsub \
-n 1 \
-q short \
-W 0:40 \
-o ${LOGS}/cleanup_%J.out \
-e ${LOGS}/cleanup_%J.err \
-w "done(${jid_reselect_som[*]}) && done(${jid_reselect_germ[*]})" \
-J "cleanup[1]" ${SCRIPT_DIR}/cleanup.lsf \
${SCRIPT_DIR} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

fi

if [[ "${RESELECT}" == "Somatic" ]] || [[ "${RESELECT}" == "Both" ]]
then

jid_reselect_som=$(bsub \
-n 4 \
-q big-multi \
-W 24:00 \
-o ${LOGS}/reselect_som_%J.out \
-e ${LOGS}/reselect_som_%J.err \
-J "reselect_som[1-${TUMOR_COUNT}]" ${SCRIPT_DIR}/reselect.lsf \
${REF} ${AIMS_DIR} ${AIMS_OUT} ${SCRIPT_DIR}/TUMOR_PATIENTS.tmp \
${COMMON_SNPS} ${DBSNP} ${SOMATIC} Somatic ${AF_DIR} \
${NORM_GERMLINE_RESOURCE} ${SCRIPT_DIR} ${NO_DECOY_INTERVALS} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

fi

if [[ "${RESELECT}" == "Germline" ]] || [[ "${RESELECT}" == "Both" ]]
then

jid_reselect_germ=$(bsub \
-n 4 \
-q big-multi \
-W 24:00 \
-o ${LOGS}/reselect_germ_%J.out \
-e ${LOGS}/reselect_germ_%J.err \
-J "reselect_germ[1]" ${SCRIPT_DIR}/reselect.lsf \
${REF} ${AIMS_DIR} ${AIMS_OUT} ${SCRIPT_DIR}/TUMOR_PATIENTS.tmp \
${COMMON_SNPS} ${DBSNP} ${GERMLINE} Germline ${AF_DIR} \
${NORM_GERMLINE_RESOURCE} ${SCRIPT_DIR} ${NO_DECOY_INTERVALS} | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

fi

# run individiual scripts here. don't forget to remove job dependencies
if [[ "${ONESCRIPT}" == "True" ]]
then

jid_CADD=$(bsub \
-n 4 \
-q big-multi \
-W 670:00 \
-o CADD_dir/CADD_%J.out \
-e CADD_dir/CADD_%J.err \
-J "CADD_som[1-${TUMOR_COUNT}]" ${SCRIPT_DIR}/CADDextracter.lsf \
${SCRIPT_DIR}/TUMOR_PATIENTS.tmp | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

fi
