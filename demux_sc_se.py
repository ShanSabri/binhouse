#!/usr/bin/env python


'''
demux_sc_se.py
Shan Sabri
Dec. 29, 2018
Demultiplex single cell ChIP-seq data (single end)
I: QSEQ
O: Compressed FQ
'''


import datetime, os, editdistance, gzip
from contextlib import contextmanager
from itertools import zip_longest as zip


def demux(iDir, rDir, idxDict, outDir, edits=2, minLen=1, verbose=True):
    """
    Demultiplex Dropseq Qseq files by splitting on sample indicies.
    """
    out = {}
    for k, v in idxDict.items():
        out[k] = gzip.open(os.path.join(outDir, v + ".fq.gz"), mode="w")

    with gzip.open(iDir, mode='r') as idxf, gzip.open(rDir, mode='r') as rf:
        for lineNum, (i, r) in enumerate(zip(idxf, rf)):

            i = i.decode('utf8').strip().split("\t")
            r = r.decode('utf8').strip().split("\t")
            obsRead = r[8].rstrip(".")
            if len(obsRead) < minLen: continue
            obsRead = obsRead.replace(".", "N")
            obsIdx = i[8]
            obsReadQual = r[9][:len(obsRead)]
            obsReadId = "@{} len:{} adp:{}".format(":".join(r[0:8]), str(len(obsRead)), obsIdx)

            def match_keys(_trueIdx):
                return editdistance.eval(_trueIdx, obsIdx)

            trueIdx = min(out.keys(), key=match_keys)

            if editdistance.eval(trueIdx, obsIdx) >= edits:
                lineNum += 1
                continue

            if verbose:
                print("\tProcessing line {}: {} {} to {}".format(lineNum, ":".join(r[0:8]), obsIdx, trueIdx))

            read = out[trueIdx]
            read.write("{0}\n{1}\n+\n{2}\n".format(obsReadId, obsRead, obsReadQual).encode('utf-8'))


@contextmanager
def log(message):
    """
    Log a timestamp, a message, and the elapsed time to console.
    """
    start = datetime.datetime.now()
    print("{} # {}".format(datetime.datetime.now(), message))
    yield
    elapsed = datetime.datetime.now() - start
    print("{} # done in {}m {}s".format(datetime.datetime.now(), divmod(elapsed.total_seconds(), 60)[0], divmod(elapsed.total_seconds(), 60)[1]))


def main():
    w = "/Users/shansabri/Desktop/"
    rDir = w + "1.qseq.txt.gz"
    iDir = w + "2.qseq.txt.gz"
    idxDict = {"ATGTCAG": "samp1"}
    outDir = w
    with log("Demultiplexing"): demux(iDir, rDir, idxDict, outDir, edits=2, minLen=1, verbose=True)


if __name__ == "__main__":
    main()