#!/bin/bash


awk '{print "chr"$1"\t"$4"\t"$4"\t"$2}' $1 > $2
