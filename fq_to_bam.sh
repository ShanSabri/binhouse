#!/bin/bash

source ~/.bash_profile

#$ -N bcs_to_bam
#$ -cwd
#$ -o /u/home/s/ssabri/project-ernst/dropseq/DROP12-map-pMX/logs
#$ -e /u/home/s/ssabri/project-ernst/dropseq/DROP12-map-pMX/logs
#$ -m bea
#$ -j y
#$ -l highp,h_data=19G,h_rt=6:00:00
#$ -t 1-4:1

echo “Task id is $SGE_TASK_ID”

java -Xmx16g -jar /u/home/s/ssabri/bin/picard-tools-1.138/picard.jar FastqToSam \
        FASTQ=/u/home/s/ssabri/project-ernst/dropseq/2018-05-16-TC-JL/demux/DROP12_${SGE_TASK_ID}_1.fq.gz \
        OUTPUT=/u/home/s/ssabri/project-ernst/dropseq/DROP12-map-pMX/DROP12_${SGE_TASK_ID}_bcs.bam \
        SAMPLE_NAME=DROP12_${SGE_TASK_ID}