# coding: utf-8

'''
demux_scrna_qseq.py
Shan Sabri
Nov 10, 2018

Demultiplex Dropseq Qseq files
I: gzip'd qseqs
O: gzip'd fqs
'''

import editdistance, gzip, os
from itertools import izip


def demux(bDir, iDir, rDir, idxDict, outDir, edits=2, minLen=1, verbose=True):
    """
    Demultiplex scRNA QSEQ files by splitting on sample indicies.

    """
    out = {}
    for k,v in idxDict.items():
        out[k] = (gzip.open(os.path.join(outDir, v + "_1.fq.gz"), "w"),
                  gzip.open(os.path.join(outDir, v + "_2.fq.gz"), "w"))

    with gzip.open(bDir) as bf, gzip.open(iDir) as idxf, gzip.open(rDir) as rf:
        for lineNum, (b, i, r) in enumerate(izip(bf, idxf, rf)):
            b = b.strip().split("\t")
            i = i.strip().split("\t")
            r = r.strip().split("\t")

            obsRead = r[8].rstrip(".")
            if len(obsRead) < minLen: continue
            obsRead = obsRead.replace(".", "N")
            obsIdx = i[8]
            obsReadQual = r[9][:len(obsRead)]
            obsBc = b[8][0:20].replace(".", "N")
            obsBcQual = b[9][0:20]

            obsReadId = "@{} len:{} bc:{}".format(":".join(b[0:8]), str(len(obsRead)), obsBc)
            obsBcId   = "@" + ":".join(b[0:8]) + " length:" + str(len(obsBc))

            def match_keys(_trueIdx):
                return editdistance.eval(_trueIdx, obsIdx)

            trueIdx = min(out.keys(), key=match_keys)

            if editdistance.eval(trueIdx, obsIdx) >= edits:
                lineNum += 1
                continue

            if verbose:
                print(
                "\tProcessing line {}: {} {} to {}".format(lineNum, ":".join(b[0:8]), obsIdx, trueIdx))

            barcode, read = out[trueIdx]
            barcode.write("{0}\n{1}\n+\n{2}\n".format(obsBcId, obsBc, obsBcQual))
            read.write("{0}\n{1}\n+\n{2}\n".format(obsReadId, obsRead, obsReadQual))