#!/bin/bash

fh=$1

head="$(head -n 1 $fh | awk '{print $1}') 1"
echo $head
grep '_' $fh
echo '#1'
echo ';'
