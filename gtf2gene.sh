#!/bin/bash
# usage: gtf2gene.sh gene.gtf ref.fa output.fa

if [ $# -ne 3 ]; then
    echo "gtf2gene.sh"
    echo "      generate combined transcripts for each gene"
    echo "Usage:"
    echo "      $0 gene.gtf ref.fa output.fa"
    exit
fi

gene_gtf=$1
ref_fa=$2
out_fa=$3
tmp=.tmp

gffread=gfftrans # modified gffread
fxtools=fxtools
pseudo_trans=./pseudo_trans.sh

# use gffread to generate all transcripts (with read name) for each gene
echo "$gffread $gene_gtf -g $ref_fa -w $tmp"
$gffread $gene_gtf -g $ref_fa -w $tmp

# use fxtools to merge transcripts of same gene
# 'N' is used to separate transcripts of one gene
echo "$fxtools merge-fa $tmp $out_fa N"
$fxtools merge-fa $tmp $out_fa N
rm $tmp

# generate pseudo longest anno-transcript for each gene, output to cluster file
bash merge_trans.sh $ref_fa $gene_gtf


# # generate bed file for only gene
# awk 'BEGIN{OFS="\t"} ($3=="gene"){print $1,$4-1,$5}' $gene_gtf > $gene_gtf.bed
# 
# bedtools=bedtools
# # use bedtools to generate sequence of each gene region
# $bedtools getfasta -fi $ref_fa -bed $gene_gtf.bed -fo $out_fa
# rm $gene_gtf.bed
