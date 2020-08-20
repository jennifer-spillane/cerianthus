## Mitochondrial genome things  

### Isolating the mitochondrial genome from the rest of the genome  

In Stampar et al. 2019, they say they isolated the mito genome by:
- creating a blast database from their whole genome assemblies (Illumina paired-end and mate-pair data assembled using discovar).
- using "known cnidarian mitochondrial coding sequences" as the query for their blast search
- mapping trimmed reads back to the blast hits
- mapped reads were reassembled (they used Geneious)

Then they did some validating and checking the stats of their new assemblies, but let's just tackle these first steps first.  

I've made a new directory within "cerianthid" to house these mitochondrial explorations: /mnt/lustre/macmaneslab/jlh1023/cerianthid/mito_explorations  

1. the above paper doesn't list which genes they used to pull out mitochondrial genes in the cerianthids, so I'm mostly looking in Kayal et al. 2013, where they have a huge table of these accession numbers. I grabbed the ones they had listed in their methods, but am also taking some from each group of anthozoans. Since this is just a first pass, I'm going to take 3-4 from each group and see how that goes.  

These are contained in the directory here: /mnt/lustre/macmaneslab/jlh1023/cerianthid/mito_explorations/mito_genes  

They are also listed with taxonomic info in a spreadsheet here (not on premise): /Users/jenniferlhill/Desktop/Analyses/cerianthid/mito_genes.xlsx  

There are some that I downloaded and they were fastq files instead of fasta, so I converted those with the fastx toolkit, like this:
`fastq_to_fasta -i SRR8109825.fastq -o SRR8109825.fa`  

2. Next I'll make a blast database out of my newly reassembled and polished cerianthus genome.  

It lives here: /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri_pol4.fasta

Making the blast database:  
`makeblastdb -in /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri_pol4.fasta -out ceri_db -dbtype nucl`  

3. And then I'll blast the mitochondrial genes against the cerianthus assembly.  

`blastn -db ceri_db -max_target_seqs 1 -query mito_genes/all_mito_genes.fa -outfmt '6 qseqid qlen length pident gaps evalue stitle' -evalue 1e-10 -num_threads 6 -out mito_blast.out`  

#### Blast results  

These are the contigs that this blast run matched to the mitochondrial genes from other cnidarians:  
- ctg2801  
- ctg1506  
- ctg1842  

2801 matched to a number of sequences, and 1506 matched a few, but 1842 only matched once, to "ENA|JX023264|JX023264.1", which is another cerianthid sequence.  

I pulled out these three contigs into their own separate files, but first had to "fix" the lines of the assembly so that they didn't wrap, which messes with my grepping.  
`awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri_pol4.fasta > unwrapped_ceri_pol4.fa`  

Then I can pull out the sequences with their headers with a context line after the match.  
`grep -A 1 "ctg2801" unwrapped_ceri_pol4.fa > ctg2801.fa`  

Each of these is a decent sized bit of sequence, 23-40K long.
