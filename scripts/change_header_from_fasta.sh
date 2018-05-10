#!/bin/bash

perl -pi -e 's/^(>.[^ ]+) .+/$1/g' $1
