#!/bin/bash
# sge_info.sh
# by Shan Sabri, github.com/ShanSabri
# 
# Description: Get usage information and availability for each node on SGE.

echo -e "queuename   \t   \t   \tresv   \tused   \ttot   \tavail"
qstat -f | grep compute | awk -F'[[:space:]]+|/' '{ OFS = "   \t"} { $6 = $5 - $4 - $3 ; print $1, $3, $4, $5, $6}'