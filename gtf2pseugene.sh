#!/bin/bash
# usage: gtf2pseu-gene.sh gene.gtf ref.fa output.fa

if [ $# -ne 3 ]; then
    echo "gtf2pseugene.sh"
    echo "      generate combined transcript for each gene"
    echo "      remove 1~2 transcript for multi-trans gene"
    echo "Usage:"
    echo "      $0 ref.fa gene.gtf pseu-gene.fa "
    exit
fi

ref_fa=$1
gene_gtf=$2
gene_fa=$3

out_dir=$(dirname $gene_fa)

exon_gtf=$out_dir/.tmp.exon.gtf
tmp=$out_dir/.tmp.fa
tmp2=$out_dir/.tmp2.fa

#gffread=gffread_GeneID_TransID # modified gffread
gffread_gene=gffread-only-GeneID
fxtools=fxtools
sort_fa=sort_fa.sh
pseudo_trans=./merge_trans.sh

# 1. extrack all 'exon' lines of gtf
echo "cat $gene_gtf | awk '\$3=="exon"{print}' > $exon_gtf"
cat $gene_gtf | awk '$3=="exon"{print}' > $exon_gtf

# 2. use gffread to generate all transcripts (with read name) for each gene
echo "$gffread_gene $exon_gtf -g $ref_fa -w $tmp"
$gffread_gene $exon_gtf -g $ref_fa -w $tmp
rm $exon_gtf

# 3. merge fasta lines
echo "$fxtools aq $tmp > $tmp2"
$fxtools aq $tmp > $tmp2
echo "$fxtools qa $tmp2 > $tmp"
$fxtools qa $tmp2 > $tmp

# 4. sort fasta with seq name
# sort $tmp
echo "$sort_fa $tmp > $tmp2"
$sort_fa $tmp > $tmp2
rm $tmp

# 5. use fxtools to merge transcripts of same gene
#    'N' is used to separate transcripts of one gene
#    if gene has more than 1 transcript, remove 1~2
echo "$fxtools merge-filter-fa $tmp2 $gene_fa N"
$fxtools merge-filter-fa $tmp2 N > $gene_fa
rm $tmp2
