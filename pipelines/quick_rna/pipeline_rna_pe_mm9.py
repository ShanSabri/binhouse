#!/usr/bin/python2
import os
from os.path import join

# CHANGE AS NEEDED 
BASE = "/u/home/k/kp1/project-kp1/Shawn_Tan/2019-08-19/"
samples = [
        "1-day-DOX-1-day-no-DOX_S5",
        "1-day-DOX-2-days-no-DOX_S8",
        "1-day-DOX_S1",
        "1-day-no-DOX_S2",
        "2-days-DOX-1-day-no-DOX_S9",
        "2-days-DOX-2-days-no-DOX_S12",
        "2-days-DOX_S3",
        "2-days-no-DOX_S4",
        "3-days-DOX-1-day-no-DOX_S13",
        "3-days-DOX-2-days-no-DOX_S16",
        "3-days-DOX_S6",
        "3-days-no-DOX_S7",
        "4-days-DOX_S10",
        "4-days-no-DOX_S11",
        "5-days-DOX_S14",
        "5-days-no-DOX_S15"
]



# GLOBALS
WRDIR   = BASE  + "process/"
LOGDIR  = WRDIR + "logs/"
DATADIR = BASE  + "data/"
GTFDIR  = "/u/home/k/kp1/project-kp1/reference_genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
REFDIR  = "/u/home/k/kp1/project-kp1/reference_genomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome"
CHRDIR  = "/u/home/k/kp1/project-kp1/reference_genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/chromInfo.txt"

SCRIPT_TEMPLATE = """\
#!/bin/bash
source ~/.bashrc
#$ -cwd
#$ -o {LOGDIR}
#$ -e {LOGDIR}
#$ -m n
#$ -l h_data=4G,h_rt=14:59:00
#$ -pe shared 6

module load python/2.7 samtools/1.9 bedtools/2.26.0 bowtie2/2.2.9

# ALIGN
~/bin/tophat2 \
        --num-threads 6 \
        --max-multihits 15 \
        --no-discordant \
        --no-mixed \
        --b2-sensitive \
        --no-novel-juncs \
	--no-novel-indels \
        --transcriptome-max-hits 1 \
        --prefilter-multihits \
        --library-type fr-firststrand \
        -G {GTFDIR} \
        -o {OUTDIR} \
        {REFDIR} \
        {MATE1} \
	{MATE2}

# SORT 
/u/local/apps/samtools/1.9/bin/samtools \
        sort \
	-@ 6 \
	-o {OUTDIR}/{S}.bam \
        {OUTDIR}/accepted_hits.bam

# INDEX
/u/local/apps/samtools/1.9/bin/samtools \
         index \
        {OUTDIR}/{S}.bam

# COUNT
~/bin/htseq-count \
        --format=bam \
        --order=pos \
        --stranded=reverse \
        --minaqual=0 \
        --type=exon \
        --mode=union \
        --idattr=gene_id \
        {OUTDIR}/{S}.bam \
        {GTFDIR} \
    > {OUTDIR}/{S}.counts

# COVERAGE 
bamCoverage \
       -p 6 \
	-of bigwig \
	--normalizeUsing=RPKM \
	-v \
	-b {OUTDIR}/{S}.bam \
	-o {OUTDIR}/{S}_rpkm.bw

# CLEAN UP 
rm {OUTDIR}/accepted_hits.bam
"""

if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)

for index, S in enumerate(samples, start=1):
    print "Processing: %s (%s/%s)" % (S, index, len(samples))

    OUTDIR = join(WRDIR, S)
    if not os.path.exists(OUTDIR): os.makedirs(OUTDIR)

    MATE1 = join(DATADIR, "%s_L001_R1_001.fastq.gz" % S)
    MATE2 = join(DATADIR, "%s_L001_R2_001.fastq.gz" % S)

    script = SCRIPT_TEMPLATE.format(**globals())

    script_name = "ST_%s" % S
    with open(script_name, 'w') as f:
        f.write(script)
    os.system("qsub %s" % script_name)
    os.system("rm %s" % script_name)

