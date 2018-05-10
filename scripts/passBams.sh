#!/bin/bash

job=$1

scp root@156.35.56.127:/datos/dreamgenics/dreampipe/jobs/$job/variantCalling/align/wrapped/postAlign/removeDuplicates/.executions/*/*/*.bam .
scp root@156.35.56.127:/datos/dreamgenics/dreampipe/jobs/$job/variantCalling/align/wrapped/postAlign/removeDuplicates/.executions/*/*/*.bam.md5 .
scp root@156.35.56.127:/datos/dreamgenics/dreampipe/jobs/$job/variantCalling/align/wrapped/postAlign/removeDuplicates/.executions/*/*/*.bam.bai .
