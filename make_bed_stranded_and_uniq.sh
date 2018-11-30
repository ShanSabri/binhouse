#!/bin/bash
# make_bed_stranded_and_uniq.sh
# by Shan Sabri, github.com/ShanSabri
# 
# Description: Take bed files from Dropseq pipeline and remove untagged reads and uniq in parallel

source ~/.bash_profile

BED_DIR=""

for f in ${BED_DIR}/*.bed.gz ; do 
	echo $f 
	gzcat $f | 
	pv | 
	parallel --pipe --keep-order --block 700M -j6 -q perl -ne 'use Digest::MD5 qw(md5_base64); print unless $seen{md5_base64($_)}++' | 
	awk '$5 == $11 {print $0}' | 
	pigz -p 6 > ${f%.*.*}.uniq.stranded.bed.gz 
done