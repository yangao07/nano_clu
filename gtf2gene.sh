#!/bin/bash
# usage: gtf2gene.sh gene.gtf ref.fa output.fa

if [ $# -ne 3 ]; then
    echo "gtf2gene.sh"
    echo "      generate combined transcripts for each gene"
    echo "Usage:"
    echo "      $0 ref.fa gene.gtf out_dir"
    exit
fi

ref_fa=$1
gene_gtf=$2
out_dir=$3
tmp=.tmp
gene_fa=$out_dir/gene.fa

gffread=gffread # modified gffread
fxtools=fxtools
pseudo_trans=./merge_trans.sh

# use gffread to generate all transcripts (with read name) for each gene
echo "$gffread $gene_gtf -g $ref_fa -w $tmp"
$gffread $gene_gtf -g $ref_fa -w $tmp

# use fxtools to merge transcripts of same gene
# 'N' is used to separate transcripts of one gene
echo "$fxtools merge-fa $tmp $gene_fa N"
$fxtools merge-fa $tmp N > $gene_fa
rm $tmp

# generate pseudo longest anno-transcript for each gene, output to cluster file
echo "bash $pseudo_trans $ref_fa $gene_gtf $out_dir"
bash $pseudo_trans $ref_fa $gene_gtf $out_dir
