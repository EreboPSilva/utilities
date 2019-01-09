#!/bin/bash

for i in $(ls *.gz);
do
        md5gz=$(md5 $i | cut -d " " -f1)
        md5su=$(cat $i.md5 | cut -d " " -f1)

        if  [ "$md5gz" = "$md5su" ]
        then
                echo "$i IS OK!"
        else
                echo "$i IS INCOMPLETE OR CORRUPTED!!! CHECK!!!"
        fi
done

