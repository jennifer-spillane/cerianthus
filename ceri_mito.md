## Mitochondrial genome things  

### Isolating the mitochondrial genome from the rest of the genome  

In Stampar et al. 2019, they say they isolated the mito genome by:
- creating a blast database from their whole genome assemblies (Illumina paired-end and mate-pair data assembled using discovar).
- using "known cnidarian mitochondrial coding sequences" as the query for their blast search
- mapping trimmed reads back to the blast hits
- mapped reads were reassembled (they used Geneious)

Then they did some validating and checking the stats of their new assemblies, but let's just tackle these first steps first.  


My genome is here: /mnt/lustre/macmaneslab/jlh1023/cerianthid/reassembly/ceri_pol7.fasta  
I've made a new directory within "reassembly" to house these mitochondrial explorations: /mnt/lustre/macmaneslab/jlh1023/cerianthid/reassembly/mito_things/

Making the blast database:  
> makeblastdb -in ceri_pol7.fasta -out mito_things/ceri_db -dbtype nucl  

That was the easy part. Now I need to find the sequences that folks expect to be on a cnidarian mitochondrial genome.  
