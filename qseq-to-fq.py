#!/usr/bin/python


__author__  = "Shan Sabri"
__email__   = "ShanASabri@gmail.com"
__date__    = "2019/03/28"


import datetime, gzip, argparse, glob, os
from itertools import izip


def demux(rd1, rd2, sample_id, output_dir):
	start_time = datetime.datetime.now()
	outfile_1 = os.path.join(output_dir, sample_id + "_1.fq.gz")
	outfile_2 = os.path.join(output_dir, sample_id + "_2.fq.gz")

	with gzip.open(rd1) as f1, gzip.open(rd2) as f2, gzip.open(outfile_1, "w") as out_1, gzip.open(outfile_2, "w") as out_2:
		print "START: ", start_time
		for line_num, (r1, r2) in enumerate(izip(f1, f2), start = 1):
			# READ1
			r1 = r1.strip().split("\t")
			r1_read = r1[8].replace(".", "N")
			r1_qual = r1[9]
			r1_id = "@" + ":".join(r1[0:8]) + " length:" + str(len(r1_read))		
			# READ2
			r2 = r2.strip().split("\t")
			r2_read = r2[8].replace(".", "N")
			r2_qual = r2[9]
			r2_id = "@" + ":".join(r2[0:8]) + " length:" + str(len(r2_read))
                	# PROGRESS
			line_num+=1
			if line_num % 10000 == 0:
				print datetime.datetime.now()-start_time, "\tProcessing line", line_num, ": {" + r1_id + ", " + r2_id + "}"
			# OUTPUT
			out_1.write("{0}\n{1}\n+\n{2}\n".format(r1_id, r1_read, r1_qual))
			out_2.write("{0}\n{1}\n+\n{2}\n".format(r2_id, r2_read, r2_qual))

	print "FINISH: ", datetime.datetime.now()


def purge(dir, pattern):
	# CLEAN UP
	for f in glob.glob(dir + pattern):
		print 'Deleting: ' + f 
		os.remove(f)


if __name__ == "__main__": 

    parser = argparse.ArgumentParser()
    parser.add_argument('read_1_file', action = "store", help = "Concatenated and compressed QSEQ file corresponding to Read 1 (e.g. read1.qseq.txt.gz)")
    parser.add_argument('read_2_file', action = "store", help = "Concatenated and compressed QSEQ file corresponding to Read 2 (e.g. read2.qseq.txt.gz)")
    parser.add_argument('sample_name', action = "store", help = "Sample name")
    parser.add_argument('output', action = "store", help = "Output directory for fastq file")
    args = parser.parse_args()

    rd1 = args.read_1_file 
    rd2 = args.read_2_file
    demux(rd1, rd2, args.sample_name, args.output)
    # purge(wrdir, "/*.qseq.txt.gz")
