#!/bin/bash

for file in G*; do
  new=$(head -1 $file | sed 's/>//');
  name=".seq";
  newname=$new$name;
  if [ -f "$newname" ];
  then
    echo "Cannot rename $file to $newname - file already exists";
  else
    mv "$file" "$newname";
  fi
done
