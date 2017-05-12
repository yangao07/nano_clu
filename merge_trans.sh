#!/bin/bash
if [[ $# -ne 2 ]]; then
    echo ""
    echo "$0 ref.fa gene.gtf"
    echo ""
    exit
fi
#ref=~/program/data/chr1.fa
ref=$1
#gtf=~/program/data/test.gtf
gtf=$2

gffread=~/software/gffread-0.9.8c/gffread
# gene.name
name=.tmp_name
cat $gtf \
    | awk -F "\"|\t| " '($3=="gene"){print $11}' > $name
# generate gene.exon.gtf
exon_bed=.tmp_exon.bed
merge_bed=.tmp_merge.bed
merge_gtf=.tmp_merge.gtf
while IFS='' read -r line || [[ -n "$line" ]]; do
    gene=$line
    # gene.exon.gtf => gene.exon.bed
    grep $gene $gtf | grep exon | awk 'BEGIN{OFS="\t"} {print $1,$4,$5}' | sort -k2,3 -n > $exon_bed
    # gene.exon.bed => merge.bed
    bedtools merge -i $exon_bed > $merge_bed
    rm $exon_bed
    # merge.bed => merge.gtf
    echo $gene
    awk -v g=$gene 'BEGIN{OFS="\t"} {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; transcript_name \"%s\";\n", $1, "NULL", "EXON", $2,$3, ".", "+", ".", g, g, g, g}' $merge_bed  > $merge_gtf
    #rm $merge_bed
    # gff merge.gtf => gene_trans.fa
    $gffread $merge_gtf -g $ref -w $gene.fa
    #rm $merge_gtf
done < "$name"
rm $name
