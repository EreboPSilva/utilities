#!/bin/bash

nfa="NM_$1.nfa";
cp $nfa seqs/.;

fa="NP_$1.fa";
transeq $nfa seqs/$fa;

rm $nfa;