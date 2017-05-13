1. Build gene annotation index
    1. input <- ref.fa, gene.gtf
    2. gene_trans.fa: generate each transcript and combined into 1 pseudo chromosome, separated by 'N'
    3. construct longest pseudo anno-transcript for each gene, output to each cluster file
    4. debwt index gene_trans.fa

2. MEM seeding and clustering
    1. count MEM seeds
    2. find most likely and secdonary likely gene
    3.  if (sec_score and best_score NOT close enougy
            1. assign read to gene
            2. output read to tmp clu file
            
        else
            assign read to remain-cluster

3. POA each cluster and corresponding pseudo anno-transcript
    1. assign anno-transcript with highest weight
    2. set gap-penalty

4. For each called consensus transcript:
    1. compare each un-aasigned read to consensus
    2. assign un-assigned read with most likely consensus/gene

5. Re-run POA for new updated clusters

6. Map consensus to the genome
    1. generate consensus' coordinate information
    2. generate all reads' coordinate information

7. Use BAM2GTF to generate updated GTF annotaion

8. Run rMATS with updated GTF to discover novel alternative splice event
