## Cerianthid assembly 2020

It was a little while ago that I did the initial assembly for this genome (before I had to set it aside for a while), so I am going to reassemble it now with updated methdods, and see if it's any better than my old one. Should be fun!  

I made a new directory to do this work in, and it is here: /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/

### The Reads  

Nanopore reads: /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/ceri234.fastq.gz  
Illumina reads: /net/storage03/backup/archive/macmanes/reads/cerianthus/  

### The Assembly  

I'll use wtdbg2 (which I guess is called Red Bean?) to do this assembly, which happens in two steps.  
Cite it here: Ruan, J. and Li, H. (2019) Fast and accurate long-read assembly with wtdbg2. Nat Methods doi:10.1038/s41592-019-0669-3  
Github here: https://github.com/ruanjue/wtdbg2   

First, the normal assembly process (using "fuzzy" de Bruijn graphs):
> /mnt/lustre/macmaneslab/macmanes/wtdbg2/wtdbg2 -x ont -g 544m -t 24 -i ceri234.fastq.gz -fo ceri_assembly

Then, the consenser:  
> /mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns -t 24 -i ceri_assembly.ctg.lay.gz -fo ceri_assembly.ctg.fa
