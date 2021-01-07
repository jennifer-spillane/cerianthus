## Cerianthid assembly 2020

It was a little while ago that I did the initial assembly for this genome (before I had to set it aside for a while), so I am going to reassemble it now with updated methods, and see if it's any better than my old one. Should be fun!  

I made a new directory to do this work in, and it is here: `/mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/`


### The Reads  

Nanopore reads: `/net/storage03/backup/archive/macmanes/reads/cerianthus/ceri234.fastq.gz`  
Illumina reads: `/net/storage03/backup/archive/macmanes/reads/cerianthus/`  

The Illumina reads are a little confusing, but I've figured out what they all are.  
The ones labeled "body" and "tentacle" are RNA-seq reads, and these read sets are both contained in the files named "ceri_rna_reads". The reads named "ceri" and "ceri_reads" are dna reads and are both subsets of the read sets named "ceri_all_reads".

I'm going to cat the DNA and RNA reads together for the polishing, and these will be called "absolutely_everything".


### The Assembly  

I'll use wtdbg2 (which I guess is called Red Bean?) to do this assembly, which happens in two steps.  
Version: 2.5 (20190621)  
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

Quast reference: QUAST: Quality assessment tool for genome assemblies  
Quast version: 4.6.0, 22f3f69   

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


**BUSCO results for some of the polishing iterations**  

BUSCO version: 3.0.0

First (ceri_pol1.fasta):  
C:66.9% [S:66.6%, D:0.3%], F:4.2%, M:28.9%, n:954

Fourth (ceri_pol4.fasta):  
C:67.2% [S:66.9%, D:0.3%], F:4.7%, M:28.1%, n:954

Fifth (ceri_pol5.fasta), run with the Metazoa dataset:  
C:87.5% [S:86.5%, D:1.0%], F:3.1%, M:9.4%, n:978

Sixth (Pachycerianthus_borealis.fa):  
Euk - C:83.5%[S:82.5%,D:1.0%],F:6.3%,M:10.2%,n:303  
Met - C:87.6%[S:86.5%,D:1.1%],F:3.0%,M:9.4%,n:978  

*gorgeous*  


I'm also going to run a python script that Dr. Adam Stuckert gave me that will give me some stats about the reads themselves. I only sequenced using nanopore for this project, despite using some other reads, so I'll just do it on those.  

First, I have to convert my reads to fasta format:  
`seqtk seq -a /net/storage03/backup/archive/macmanes/reads/cerianthus/ceri234.fastq.gz > ceri_nano_reads.fa`  

Then I can run the script that Adam wrote on my new fasta file:  
`/mnt/lustre/macmaneslab/shared/scripts/readlengths.py ceri_nano_reads.fa ceri_nano_read_lengths.txt`  

And the script outputs a textfile that has all the lengths of all the reads in it, but also pops out to the screen some stats:  
Total reads in file: 3572319  
Total bases in file: 17660674857  
Average read length: 4943.756382618685  
Read N50: 7682  

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

3. I'll need a list of all the contig names in the genome, and they can't have the normal ">" in front of them, so I'll remove that. I'll also remove the length info that is sitting behind the names.  
`grep ">" ceri_assembly.ctg.fa | sed 's_>__' | awk '{print $1}' > contig_names1.txt`  

4. Then I can split the contig names into chunks that I can run on multiple machines.  
`shuf contig_names1.txt | split -d -l 200 - genomechunk.`

5. The first 10 genome chunks will have to be renamed because pilon doesn't like the double digit number. I'll also take this opportunity to move them into a new directory called "chunks1" to keep things more organized.  
`mkdir chunks1
rename genomechunk.0 genomechunk. genomechunk.0*
mv genomechunk* chunks1/`

6. Originally I had steps 3-5 in the script with pilon also, but this ends up being problematic because of the way that the pilon job gets split up and run again and again. So I just do those steps manually, and then the pilon step is in a script called "pilon1.sh".   

Pilon docs: https://github.com/broadinstitute/pilon/wiki  
Pilon citation: Bruce J. Walker, Thomas Abeel, Terrance Shea, Margaret Priest, Amr Abouelliel, Sharadha Sakthikumar, Christina A. Cuomo, Qiandong Zeng, Jennifer Wortman, Sarah K. Young, Ashlee M. Earl (2014) Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. PLoS ONE 9(11): e112963. doi:10.1371/journal.pone.0112963  

And now I can do the actual polishing step (pilon1.sh).   
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


7. Pilon will spit out it's own log files for every chunk, as well as output files that start with "pilonchunk" and have numbers that correspond with the genomechunk numbers. These should be catted together at the end, to form the next version of the genome:  
`cat pilonchunk.* > ceri_pol1.fasta`  

8. Then I'm also going to move all the individual pilonchunk output files into a new directory so that they aren't just everywhere.  
`mkdir pilon_output1
mv pilonchunk.* pilon_output1/`  

And now the process can start over again.
Starting up at the "Polishing" heading, with the new version of the assembly.
Should be run until BUSCO values start to level off and we aren't really gaining any new info.

**annoying thing to watch out for**   
At one point during the polishing stage, there was something wrong with the genomechunks I made. A couple of them had headers inside that had strings of nucleotides in front of them, as if the previous sequence had not had a linefeed at the end, and when I grepped for the > symbol, it was grabbing these sequences also. Then pilon would choke on that entire genomechunk file, and give me back only a subset of results files, which obviously was not going to work. I tried removing those strings of nucleotides, but when I reran pilon, the full assembly still wasn't the size that it should have been, even though there were the correct number of output files. So some sequences were getting lost somewhere.  

In the end, I made new genomechunks in exactly the same way as before, and did not have the same problem again. This is not a very satisfactory conclusion, but it is working again, so hopefully if the same thing happens later in the polishing phase, I can just do the same thing.  

**8-21-20 update**  
I am going to polish once more. I thought that I needed both the forward and reverse reads to polish correctly, but apparently this is not the case, and the RNA-seq data for this organism is heavily biased in favor of forward reads. So I'm going to add those in to the mapping step, generate a new bam file (ceri5_sorted.bam). Just to make sure I'm not missing anything, I'm going to add the total DNA reads (ceri_all_reads) and these new RNA-seq read files (CB_all_rna) and pop them in files called "all_dna_and_rna".  

**11-13-20 update**  
Going to polish again (obviously). A paper came out that had some cerianthid transcriptomes in it, including our species. I tried to use these just for the maker run to get a better annotation, but we're still not where we're hoping to be with busco scores, so I'm going to try to use the RNAseq reads from our species in a polishing step (even though it's not the same individual, it seems like a solid bet that it will make it a bit better anyway). It's definitely worth a shot, even if it doesn't end up amounting to a whole lot of difference. So I added the reads onto the bottom of the "all_dna_and_rna" files (`cat Pachycerianthus_borealis_1.fastq >> /net/storage03/backup/archive/macmanes/reads/cerianthus/all_dna_and_rna1.fq.gz` for both of them), and did another bwa run and pilon run. Busco scores can be found above with the others when it finishes.  


#### The current (11-14-20) "final" iteration of this assembly is here: /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/Pachycerianthus_borealis.fa  


### Mapping  

BUSCO is in the process of being updated, so those results will be forthcoming. In the meantime, I'll look at mapping rates to try to see the polishing quality level off a bit.  

All of these get run with a really basic command, just replacing the argument with the bam file of choice.  

Samtools version: 1.10  
Samtools ref:

**Mapping the Illumina reads to the raw assembly (ceri_assembly.ctg.fa)**  
`samtools flagstat ceri1.sorted.bam`  
416799478 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
54998830 + 0 supplementary
0 + 0 duplicates
413696158 + 0 mapped (99.26% : N/A)
361800648 + 0 paired in sequencing
180900324 + 0 read1
180900324 + 0 read2
320358916 + 0 properly paired (88.55% : N/A)
357486390 + 0 with itself and mate mapped
1210938 + 0 singletons (0.33% : N/A)
34328800 + 0 with mate mapped to a different chr
18734079 + 0 with mate mapped to a different chr (mapQ>=5)  

**Mapping the Illumina reads to the first polishing iteration (ceri_pol1.fasta)**  
`samtools flagstat ceri2.sorted.bam`  
413531825 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
51731177 + 0 supplementary
0 + 0 duplicates
410750464 + 0 mapped (99.33% : N/A)
361800648 + 0 paired in sequencing
180900324 + 0 read1
180900324 + 0 read2
325200186 + 0 properly paired (89.88% : N/A)
358001758 + 0 with itself and mate mapped
1017529 + 0 singletons (0.28% : N/A)
30026094 + 0 with mate mapped to a different chr
16015076 + 0 with mate mapped to a different chr (mapQ>=5)  

**Mapping the Illumina reads to the second polishing iteration (ceri_pol2.fasta)**  
`samtools flagstat ceri3.sorted.bam`  
412327986 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
50527338 + 0 supplementary
0 + 0 duplicates
409567849 + 0 mapped (99.33% : N/A)
361800648 + 0 paired in sequencing
180900324 + 0 read1
180900324 + 0 read2
325995352 + 0 properly paired (90.10% : N/A)
358037952 + 0 with itself and mate mapped
1002559 + 0 singletons (0.28% : N/A)
29243086 + 0 with mate mapped to a different chr
15742789 + 0 with mate mapped to a different chr (mapQ>=5)

**Mapping the Illumina reads to the third polishing iteration**  
`samtools flagstat ceri4.sorted.bam`  
412031446 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
50230798 + 0 supplementary
0 + 0 duplicates
409275511 + 0 mapped (99.33% : N/A)
361800648 + 0 paired in sequencing
180900324 + 0 read1
180900324 + 0 read2
326283200 + 0 properly paired (90.18% : N/A)
358045072 + 0 with itself and mate mapped
999641 + 0 singletons (0.28% : N/A)
28958524 + 0 with mate mapped to a different chr
15722413 + 0 with mate mapped to a different chr (mapQ>=5)

**Mapping the Illumina reads to the fourth polishing iteration - this one also has extra RNA reads**  
`samtools flagstat ceri4.sorted.bam`  
411402714 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
49968231 + 0 supplementary
0 + 0 duplicates
408651647 + 0 mapped (99.33% : N/A)
361434483 + 0 paired in sequencing
180717242 + 0 read1
180717241 + 0 read2
326128285 + 0 properly paired (90.23% : N/A)
357686407 + 0 with itself and mate mapped
997009 + 0 singletons (0.28% : N/A)
28754970 + 0 with mate mapped to a different chr
15655026 + 0 with mate mapped to a different chr (mapQ>=5)
