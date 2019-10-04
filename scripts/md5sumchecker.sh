#!/bin/bash

if test -z "$1"
    then
        echo "Please specify the container format, e.g. gz, zip, b2, ..."
        exit 1
fi

for i in $(ls *.$1);
do
        md5gz=$(md5sum $i | cut -d " " -f1)
        md5su=$(cat $i.md5 | cut -d " " -f1)

        if  [ "$md5gz" = "$md5su" ]
        then
                echo "$i IS OK!"
        else
                echo "$i IS INCOMPLETE OR CORRUPTED!!! CHECK!!!"
        fi
done

