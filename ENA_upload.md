## Uploading reads and assemblies to the ENA  

I need to deposit all of my Pachycerianthus reads and assemblies and everything into the ENA to make them publicly available ahead of the genome note submission. When I inevitably need to do this again, these are my notes of what I did.  

### Preparing the files  

These are the files I need to deposit:  
- Illumina read files: /net/storage03/backup/archive/macmanes/reads/cerianthus/all_dna_and_rna1.fq.gz  
                       /net/storage03/backup/archive/macmanes/reads/cerianthus/all_dna_and_rna2.fq.gz  
- Nanopore read file: /net/storage03/backup/archive/macmanes/reads/cerianthus/ceri234.fastq.gz  
- Polished assembly: /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/Pachycerianthus_borealis.fa  
- Protein predictions: /mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/Pachycerianthus_borealis_proteins.fa  
- Gff file: /mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/pachyceri_maker.gff3  

All files should be zipped ahead of time. All my read files were, but I zipped the assembly and other files.  
`gzip /mnt/lustre/macmaneslab/jlh1023/cerianthid/assembly_2020/Pachycerianthus_borealis.fa`  

*Update: I realize I actually need the nanopore fast5 files, which I don't have access to currently. Hopefully I can figure out a way to make that happen soon. I also pushed the Illumina data onto Dave's plate, because he has the original files and knows more about how they were sampled and generated.*  

Focusing on the genome for now, since I can do that part. I need to combine the genome and annotation into a single file (called a flat file), and I can do that with a program called EMBLmyGFF3. It took me forever to figure out how to install this program, but I finally got the yml file from Adam that he used, and I could do it with that. I moved all the relevant files into a new directory for this, as it seemed easier at the time. It is here: /mnt/lustre/macmaneslab/jlh1023/embl  
    `conda env create -f myEMBLmyGFF3.yml`  
    `source activate EMBLmyGFF3`  

Then I put all of the different options into a script called `options.sh` because there are quite a few of them, and it seemed easier to be able to run it again later. It took a few hours, and finished the conversion with the file `Pachycerianthus_borealis_genomeassembly_1.0.embl`.  

*Another update, because Dave started uploading to the SRA, and not the ENA (he claims it is easier), I now have to do that instead. Unfortunately, they have a completely different format for annotated genomes, so I'll have to figure out a different way to do this.*  
