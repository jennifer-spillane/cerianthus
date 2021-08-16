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

2801 matched to a number of sequences, and 1506 matched a few, but 1842 only matched once, to "ENA|JX023264|JX023264.1 Ceriantheopsis americana NADH dehydrogenase subunit 4 (nad4) gene, partial cds; mitochondrial.".  

I pulled out these three contigs into their own separate files, but first had to "fix" the lines of the assembly so that they didn't wrap, which messes with my grepping.  
`awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri_pol4.fasta > unwrapped_ceri_pol4.fa`  

Then I can pull out the sequences with their headers with a context line after the match.  
`grep -A 1 "ctg2801" unwrapped_ceri_pol4.fa > ctg2801.fa`  

Each of these is a decent sized bit of sequence, 23-40K long.  

4. Trying blast again but this time with all individual mitochondrial genes, instead of whole mito genomes. Just want to check if I can get any different scaffolds.  

I remade the blast database, because I have a slightly updated polishing iteration of the genome:  
`makeblastdb -in /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri_pol5.fasta -out ceri_db2 -dbtype nucl`  

Then blasted all the individual mito genes against it:  
`blastn -db ceri_db2 -max_target_seqs 1 -query ind_mito_genes/all_mito_genes.fa -outfmt '6 qseqid qlen length pident gaps evalue stitle' -evalue 1e-10 -num_threads 6 -out mito_blast2.out`  

**I got the exact same contigs when I blasted this way, so I'm going to say that this is what we've got unless we think of something else to try**

#### Mapping Explorations  

I catted the three contigs that I pulled out from the blast results into a single file called "putative_mito_genome.fa" and then tried mapping both the Nanopore reads and Illumina reads back to it.  

These mapping rates are both terrible, but that was sort of to be expected.  

`samtools flagstat mito_nano.sorted.bam`  
3575675 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
3356 + 0 supplementary
0 + 0 duplicates
9709 + 0 mapped (0.27% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)  

`samtools flagstat mito_illumina.sorted.bam`  
361401721 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
25007 + 0 supplementary
0 + 0 duplicates
2473248 + 0 mapped (0.68% : N/A)
361376714 + 0 paired in sequencing
180688357 + 0 read1
180688357 + 0 read2
2412128 + 0 properly paired (0.67% : N/A)
2417998 + 0 with itself and mate mapped
30243 + 0 singletons (0.01% : N/A)
3496 + 0 with mate mapped to a different chr
3065 + 0 with mate mapped to a different chr (mapQ>=5)


### Looking at isolation a different way  

Since the mito genome was sequenced at the same time as the nuclear genome, but there were many more copies, it is very possible that it has vastly higher coverage than most/all of the nuclear genes. So I'm going to look at coverage as a way of potentially finding pieces of the mito genome.  

First, I calculated the depth of the coverage at every position in the genome assembly with this command:  
`samtools depth ceri5.sorted.bam > depth_ceri_pol5.tsv`  

Then I realized that this files was way too hilariously large to do much with (17Gb), so I decided to pull out just certain parts of it. I rediscovered a handy script I must have written ages ago that pulls out all of the contigs in a fasta file that are above a certain length. I used this to retain contigs over 1000bp.  
`/mnt/lustre/macmaneslab/jlh1023/pipeline_dev/pipeline_scripts/long_fasta.py -i ceri_pol5.fasta -o long_contigs_ceripol5.fa -l 1000`  

Unfortunately (or fortunately?) it's not clear that there are any contigs that got excluded based on this length threshold, as the resulting file is the same size as the original polished assembly from which it was made. So I subsampled down to 1% of the original file size by only grabbing every 100th line of the file.  
`/mnt/lustre/macmaneslab/jlh1023/pipeline_dev/pipeline_scripts/sub_file.py -i depth_ceri_pol5.tsv -o sub_depth.tsv -n 100`  

This still results in a file that is 165Mb, but it's still a little more manageable. I downloaded this file and changed it to a csv in bbedit (did not trust excel with a file this size) and gave it a header line, and then popped it into R for plotting. The resulting density plot was not really that satisfying (/Users/jenniferlhill/Desktop/Analyses/cerianthid/seq_depth_density.pdf). I did it first with all of the lines of data, but it went out so far as to make the rest unreadable (~8000). So I limited the x axis to 2000, and it was a little better, and then I limited it to 1000 so I could actually see what was going on in the rest of it. Sadly, there doesn't seem to be a peak of any kind at the higher depth places, it just trails off gradually forever, so I'm not sure I'll be able to find the mito contigs like this. Still, it might be worth quantifying a little what is going on with the depth at those higher points.  

I wrote a script that pulls out depth info for entries at or above a user-chosen depth threshold, and ran it:  
`/mnt/lustre/macmaneslab/jlh1023/pipeline_dev/pipeline_scripts/parse_depth.py -i sub_depth.tsv -o deep_contigs.tsv -d 1000`  

There are 4,188,497 lines of information in this sub_depth.tsv file, and 32,241 have sequence depth of 1000 or more. I'm going to up it a bit more to see if we can narrow things a bit.   

Threshold = 1500; lines = 14337  
Threshold = 2000; lines = 8389  
Threshold = 2500; lines = 5941  

I know from the plot that it will go pretty high, so I'm going to try one pretty far out.  
Threshold = 8000; lines = 115  
At 8500 and 9000 there aren't any, so I'll take a closer look at the ones I have.  

Unfortunately they look like they are from a ton of different contigs, so this is probably not going to be overly helpful for mito contig isolation. I also ran it with the threshold at 8100 (the largest one I could get results for) and it returned 3 contigs, none of which matches the ones I found blasting, so I think this is a dead end.  




### Trying out GetOrganelle  

Dave recently told me about a new program called GetOrganelle, and I thought I would play with it and see if I can pull out the mito genome of our cerianthid this way.  

To use on premise:  
`module load anaconda/colsa`  
`conda find get_organelle_from_assembly.py -a`  
or
`conda activate getorganelle-1.7.4.1`  

This is the GetOrganelle example script, so I know I need at least these parts.  
`get_organelle_from_reads.py -1 forward.fq -2 reverse.fq -o plastome_output -R 15 -k 21,45,65,85,105 -F embplant_pt`  

This is the script I will try first:  
`get_organelle_from_reads.py -1 /net/storage03/backup/archive/macmanes/reads/cerianthus/ceri_reads_1.fastq.gz -2 /net/storage03/backup/archive/macmanes/reads/cerianthus/ceri_reads_2.fastq.gz -u /net/storage03/backup/archive/macmanes/reads/cerianthus/ceri234.fastq.gz -o getorganelle_output -F animal_mt`  

I immediately got an error "ERROR: default animal_mt database not added yet!", so I used their supplied command to install that database.  
`get_organelle_config.py -a animal_mt`  

When I ran it again, I got an error asking me to check the integrity of my input reads. I'm going to try it just with the illumina ones, to see if the nanopore ones are the problem.  

This time it errored because my output directory already existed. Running again. And now I'm still getting the error about my reads. I'll try just the nanopore ones and see if they are any better.  
`get_organelle_from_reads.py -u /net/storage03/backup/archive/macmanes/reads/cerianthus/ceri234.fastq.gz -o getorganelle_output -F animal_mt`  

Seems to be running! Those sketchy as Illumina reads appear to have struck again. We'll see if these nanopore ones pan out.  



