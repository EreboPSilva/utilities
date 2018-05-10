#!/bin/bash

for i in $(ls *.gz);
do
        md5gz=$(md5sum $i)
        md5su=$(cat $i.md5)

        if  [ "$md5gz" == "$md5su" ]
        then
                echo "$i IS OK!"
        else
                echo "$i IS INCOMPLETE OR CORRUPTED!!! CHECK!!!"
        fi
done

