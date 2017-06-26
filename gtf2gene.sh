#!/bin/bash
# usage: gtf2-gene.sh gene.gtf ref.fa output.fa

if [ $# -ne 4 ]; then
    echo "gtf2-gene.sh"
    echo "      generate combined transcript for each gene"
    echo "      generate pseudo transcript for each gene (combined exonic sequences)"
    echo "Usage:"
    echo "      $0 ref.fa gene.gtf gene.fa exonic.fa"
    exit
fi

ref_fa=$1
gene_gtf=$2
gene_fa=$3
exonic_fa=$4

out_dir=$(dirname $gene_fa)

exon_gtf=$out_dir/.tmp.exon.gtf
tmp=$out_dir/.tmp.fa
tmp2=$out_dir/.tmp2.fa

gffread=gffread-only-GeneID # modified gffread
fxtools=fxtools
sort_fa=sort_fa.sh
pseudo_trans=./merge_trans.sh

# extrack all 'exon' lines of gtf
echo "cat $gene_gtf | awk '\$3=="exon"{print}' > $exon_gtf"
cat $gene_gtf | awk '$3=="exon"{print}' > $exon_gtf

# use gffread to generate all transcripts (with read name) for each gene
echo "$gffread $exon_gtf -g $ref_fa -w $tmp"
$gffread $exon_gtf -g $ref_fa -w $tmp

echo "$fxtools aq $tmp > $tmp2"
$fxtools aq $tmp > $tmp2
echo "$fxtools qa $tmp2 > $tmp"
$fxtools qa $tmp2 > $tmp

# sort $tmp
echo "$sort_fa $tmp > $tmp2"
$sort_fa $tmp > $tmp2
rm $tmp

# use fxtools to merge transcripts of same gene
# 'N' is used to separate transcripts of one gene
echo "$fxtools merge-fa $tmp2 $gene_fa N"
$fxtools merge-fa $tmp2 N > $gene_fa
rm $tmp2

# generate pseudo longest anno-transcript for each gene, output to cluster file
echo "bash $pseudo_trans $ref_fa $gene_gtf $exon_gtf $exonic_fa"
bash $pseudo_trans $ref_fa $gene_gtf $exon_gtf $exonic_fa

rm $exon_gtf
