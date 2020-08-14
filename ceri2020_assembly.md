## Cerianthid assembly 2020

It was a little while ago that I did the initial assembly for this genome (before I had to set it aside for a while), so I am going to reassemble it now with updated methdods, and see if it's any better than my old one. Should be fun!  

I made a new directory to do this work in, and it is here: `/mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/`


### The Reads  

Nanopore reads: `/mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri234.fastq.gz`  
Illumina reads: `/net/storage03/backup/archive/macmanes/reads/cerianthus/`  

The Illumina reads are a little confusing, but I've figured out what they all are.  
The ones labeled "body" and "tentacle" are RNA-seq reads, and these read sets are both contained in the files named "ceri_rna_reads". The reads named "ceri" and "ceri_reads" are dna reads and are both subsets of the read sets named "ceri_all_reads".

I'm going to cat the DNA and RNA reads together for the polishing, and these will be called "absolutely_everything".


### The Assembly  

I'll use wtdbg2 (which I guess is called Red Bean?) to do this assembly, which happens in two steps.  
Cite it here: Ruan, J. and Li, H. (2019) Fast and accurate long-read assembly with wtdbg2. Nat Methods doi:10.1038/s41592-019-0669-3  
Github here: https://github.com/ruanjue/wtdbg2   

First, the normal assembly process (using "fuzzy" de Bruijn graphs):
`/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtdbg2 -x ont -g 544m -t 24 -i ceri234.fastq.gz -fo ceri_assembly`  

I'm using 544mb as the estimated genome size, since that's what I got when I assembled it before using flye.

Then, the consenser:  
`/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns -t 24 -i ceri_assembly.ctg.lay.gz -fo ceri_assembly.ctg.fa`

That finished overnight some time (started midday-ish), so now I'll move on to some quality checks!
Just for reference, the raw assembly is here: /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri_assembly.ctg.fa


### Quality Checks

First I'm running quast, as it should give me some good contiguity numbers without being polished. I'll check them again after polishing to make sure nothing changed that dramatically.

Quast command:
`quast ceri_assembly.ctg.fa -o quast_raw_assembly -t 24`

Quast results:
/mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/quast_raw_assembly/report.txt

Largest contig        2254400
Total length          491647159
GC (%)                35.44
N50                   396334  


I'm also going to run assemblathon on this assembly, which should give me a few more things.  

Assemblathon command:  
`/mnt/lustre/macmaneslab/macmanes/bin/assemblathon_stats.pl ceri_assembly.ctg.fa > assemblathon_raw_assembly.txt`

Assemblathon results:  
/mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/assemblathon_raw_assembly.txt  

Number of scaffolds       5833
Total size of scaffolds  491647159
Longest scaffold    2254400
Number of scaffolds > 100K nt       1071  18.4%
Number of scaffolds > 1M nt         48   0.8%
Mean scaffold size      84287
Median scaffold size      15663
N50 scaffold length     396334


Next I'll run busco, which I'm expecting to go pretty poorly. First, because running busco on a genome is a huge pain, and second, because these are error-prone nanopore reads, so it probably won't be able to recognize much.



### Polishing  

I'll use the Illumina reads that we have (genomic and transcriptomic) to polish the error-y nanopore reads. I'm also pretty much stealing this procedure from when I did it before, but with the updated assembly.  

To make sure I'm not taking up any more space than I need to, I'm deleting any sam files at the end immediately, and getting rid of the bam files as soon as I'm finished polishing.  

These first two steps will be run in a script called "bwa1.sh", with each iteration having a different number.

BWA version: 0.7.17-r1188  
BWA docs: http://bio-bwa.sourceforge.net/bwa.shtml  
BWA citation: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]  

1. I'm running bwa to get an index of the assembly  
`bwa index -p index1 ceri_assembly.ctg.fa`

The "index1" is what the output file will be called, and the "ceri_assembly.ctg.fa" is what I want it to index.

2. And then aligning the reads using bwa  
`bwa mem -t 24 index1 \
/net/storage03/backup/archive/macmanes/reads/cerianthus/absolutely_everything_1.fq.gz \
/net/storage03/backup/archive/macmanes/reads/cerianthus/absolutely_everything_2.fq.gz \
| samtools view -@20 -Sb - \
| samtools sort -T ceri -O bam -@20 -l9 -m2G -o ceri1.sorted.bam -
samtools index ceri1.sorted.bam`

I'll need the name of the index, "index1", the reads that I am going to use "absolutely_everything_1.fq.gz" (and the reverse ones), then the name of the output (twice) "ceri1_sorted.bam".  

3. Polishing using Pilon. The prelims for polishing and the polishing itself are contained in a script called "pilon1.sh".  

Pilon docs: https://github.com/broadinstitute/pilon/wiki  
Pilon citation: Bruce J. Walker, Thomas Abeel, Terrance Shea, Margaret Priest, Amr Abouelliel, Sharadha Sakthikumar, Christina A. Cuomo, Qiandong Zeng, Jennifer Wortman, Sarah K. Young, Ashlee M. Earl (2014) Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. PLoS ONE 9(11): e112963. doi:10.1371/journal.pone.0112963  

I'll need a list of all the contig names in the genome, and they can't have the normal ">" in front of them, so I'll remove that. I'll also remove the length info that is sitting behind the names.  
`grep ">" ceri_assembly.ctg.fa | sed 's_>__' | awk '{print $1}' > contig_names1.txt`  


Then I can split the contig names into chunks that I can run on multiple machines.  
`shuf contig_names1.txt | split -d -l 200 - genomechunk.`

The first 10 genome chunks will have to be renamed because pilon doesn't like the double digit number. I'll also take this opportunity to move them into a new directory called "chunks1" to keep things more organized.  
`mkdir chunks1
rename genomechunk.0 genomechunk. genomechunk.0*
mv genomechunk* chunks1/`

And now I can do the actual polishing step.   
`echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID`

`java -jar -Xmx100G /mnt/lustre/macmaneslab/shared/spillane/pilon-1.23.jar \
--genome /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri_assembly.ctg.fa \
--bam /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri1.sorted.bam \
--output pilonchunk.$SLURM_ARRAY_TASK_ID \
--fix bases,gaps \
--diploid \
--threads 24 \
--flank 5 \
--verbose \
--mingap 1 \
--nostrays \
--targets /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/chunks1/genomechunk.$SLURM_ARRAY_TASK_ID`

Pilon will spit out it's own log files for every chunk, as well as output files that start with "pilonchunk" and have numbers that correspond with the genomechunk numbers. These should be catted together at the end, to form the next version of the genome:  
`cat pilonchunk.* > ceri_pol1.fasta`  

Then I'm also going to move all the individual pilonchunk output files into a new directory so that they aren't just everywhere.  
`mkdir pilon_output1
mv pilonchunk.* pilon_output1/`  

And now the process can start over again.
Starting up at the "Polishing" heading, with the new version of the assembly.
Should be run until BUSCO values start to level off and we aren't really gaining any new info.
