#!/usr/bin/python
import os
from os.path import join

WRDIR   = "/u/home/s/ssabri/project-ernst/analysis/Kostas-Chronis/2019-08-05/process/"
DATADIR = "/u/home/s/ssabri/project-ernst/analysis/Kostas-Chronis/2019-08-05/data/"
LOGDIR  = WRDIR + "logs"
REFDIR  = "/u/home/s/ssabri/project-ernst/analysis/Kostas-Chronis/2019-08-05/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"

SCRIPT_TEMPLATE = """\
#!/bin/bash
source ~/.bash_profile
#$ -cwd
#$ -o {LOGDIR}
#$ -e {LOGDIR}
#$ -m n
#$ -N ss_{s}
#$ -l highp,h_data=4G,h_rt=16:00:00
#$ -pe shared 8

# ALIGN
~/bin/bowtie2 \
       -p 8 \
       --qc-filter \
       -x  {REFDIR} \
       -1 {MATE1} \
       -2 {MATE2} \
       | samtools view -bhSo {outdir}/{s}.bam

# SORT
samtools sort \
  -@ 8 \
  -m 2G \
  -o {outdir}/{s}_sorted.bam \
  {outdir}/{s}.bam

# INDEX 
samtools index \
  {outdir}/{s}_sorted.bam
  
# BW
bamCoverage \
  -p 8 \
  -of bigwig \
  --normalizeUsing=RPKM \
  -v \
  -b {outdir}/{s}_sorted.bam \
  -o {outdir}/{s}_rpkm.bw
"""

if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)

samples = [line.strip() for line in open(WRDIR + "/samples_hg38.txt")]

for index, s in enumerate(samples, start=1):
    print "Processing: %s (%s/%s)" % (s, index, len(samples))

    outdir = join(WRDIR, s)
    if not os.path.exists(outdir): os.makedirs(outdir)

    MATE1 = join(DATADIR, "%s_R1_001.fastq.gz" % s)
    MATE2 = join(DATADIR, "%s_R2_001.fastq.gz" % s)

    script = SCRIPT_TEMPLATE.format(**globals())

    script_name = "%s" % s
    with open(script_name, 'w') as f:
        f.write(script)
    os.system("qsub %s" % script_name)
    os.system("rm %s" % script_name)
