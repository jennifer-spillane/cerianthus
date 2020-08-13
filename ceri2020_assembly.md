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

1. I'm running bwa to get an index of the assembly
`bwa index -p ceri_index1 ceri_assembly.ctg.fa`

The "ceri_index" is what the output file will be called, and the "assembly_flye/scaffolds.fasta" is what you want it to index.
This is contained in the script located here:
/mnt/lustre/macmaneslab/jlh1023/cerianthid/reassembly/index.sh

2. And then aligning the reads using bwa
>bwa mem -t 24 \
ceri_index \
/mnt/lustre/macmaneslab/jlh1023/cerianthid/ceri_illreads/ceri_1.fastq.gz \
/mnt/lustre/macmaneslab/jlh1023/cerianthid/ceri_illreads/ceri_2.fastq.gz \
| samtools view -@20 -Sb - \
| samtools sort -T ceri -O bam -@20 -l9 -m2G -o ceri.sorted.bam -
samtools index ceri.sorted.bam

You'll need the name of the index, "ceri_index", the reads that you are going to use "ceri_1.fastq.gz" and the reverse ones, then the name of the output (twice) "ceri_sorted.bam"
Contained in the script located here:
/mnt/lustre/macmaneslab/jlh1023/cerianthid/reassembly/get_bam.sh


get rid of sam files immediately
need bam files for polishing

zcat ceri_all_reads_1.fastq.gz ceri_rna_reads_1.fq.gz > absolutely_everything_1.fq.gz

ceri_all_reads_1.fastq.gz <- used last time I polished.
