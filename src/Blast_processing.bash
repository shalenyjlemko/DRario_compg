#!/bin/bash

# ========================================================
# Blast processing for Zebrafish proteome results into MCL clustering
# ========================================================


BLAST="zebrafish_blast_results_clean.txt"

awk -v OFS='\t' '!/^#/{
  ql=$13; sl=$14; qs=$7; qe=$8; ss=$9; se=$10;
  if(ql>0 && sl>0){
    qc=((qe>=qs?qe-qs+1:qs-qe+1)/ql)*100;
    sc=((se>=ss?se-ss+1:ss-se+1)/sl)*100;
    if(qc>=80 && sc>=80 && $3>=70 && $11<=1e-10 && $12>=80)
      printf "%s\t%.2f\t%.2f\n",$0,qc,sc
  }}' "$BLAST" > dataset_stringent.txt

  wc -l dataset_stringent.txt

  awk -v OFS='\t' '{print $1,$2,$12}' dataset_stringent.txt > zebrafish_mcl_input.txt