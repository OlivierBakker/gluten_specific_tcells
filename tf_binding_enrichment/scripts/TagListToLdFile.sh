#!/bin/bash



tag_list=$1
out=$2

awk '{print $1":\t"$1"|"$8}' $tag_list | tail -n +2 | sed 's/|/\t/g' | sed 's/NONE//g' > $2
