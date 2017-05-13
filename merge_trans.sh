#!/bin/bash
if [[ $# -ne 3 ]]; then
    echo ""
    echo "$0 ref.fa gene.gtf out_dir"
    echo ""
    exit
fi
ref=$1
gtf=$2
out_dir=$3

gffread=gffread
# gene.name
name=$out_dir/.tmp_name
cat $gtf | awk -F "\"|\t| " '($3=="gene"){print $11}' > $name
# generate gene.exon.gtf
exon_bed=$out_dir/.tmp_exon.bed
merge_bed=$out_dir/.tmp_merge.bed
merge_gtf=$out_dir/.tmp_merge.gtf
while IFS='' read -r line || [[ -n "$line" ]]; do
    gene=$line
    out=$out_dir/$gene.fa
    # gene.exon.gtf => gene.exon.bed
    grep $gene $gtf | grep exon | awk 'BEGIN{OFS="\t"} {print $1,$4,$5}' | sort -k2,3 -n > $exon_bed
    # gene.exon.bed => merge.bed
    bedtools merge -i $exon_bed > $merge_bed
    rm $exon_bed
    # merge.bed => merge.gtf
    awk -v g=$gene 'BEGIN{OFS="\t"} {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\";\n", $1, "NULL", "EXON", $2,$3, ".", "+", ".", g, "NULL", g}' $merge_bed  > $merge_gtf
    rm $merge_bed
    # gff merge.gtf => gene_trans.fa
    echo "$gffread $merge_gtf -g $ref -w $out"
    $gffread $merge_gtf -g $ref -w $out
    rm $merge_gtf
done < "$name"
rm $name
