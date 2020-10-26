## Running MAKER on the new Cerianthus assembly

I need to run MAKER on my new assembly, so first I am gathering all of the other Cnidarians that I will use in the run. All of these will be in this directory: `/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/`  

Maker version: 3.01.02  
Maker docs: http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#MAKER.27s_Output  
Maker ref: MAKER: an easy-to-use annotation pipeline designed for emerging model organism genomes.  

### First run of MAKER  

**Proteins:**  
- Acropora_digitifera_prot.fa  
- Acropora_millepora_prot.fa  
- Hydra_magnipapillata_prot.fa  
- Nematostella_vectensis_prot.fa  
- Thelohanellus_kitauei_prot.fa  

**Transcriptomes**  
- Clytia_hemisphaerica_trans.fa  
- Hydra_vulgaris_trans.fa  
- Kudoa_iwatai_ORP_223_TMP.orthomerged.fasta  
- Alatina_alata_ORP_223_TMP.orthomerged.fasta  
- Liriope_tetraphylla_ORP_223_TMP.orthomerged.fasta  
- Alcyonium_palmatum_ORP_223_TMP.orthomerged.fasta  
- Lucernariopsis_campanulata_ORP_223_TMP.orthomerged.fasta  
- Antipathes_caribbeana_ORP_223_TMP.orthomerged.fasta  
- Myxobolus_cerebralis_ORP_223_TMP.orthomerged.fasta  
- Chironex_fleckeri_ORP_223_TMP.orthomerged.fasta  
- Pelagia_noctiluca_ORP_223_TMP.orthomerged.fasta  
- Corallium_rubrum_ORP_223_TMP.orthomerged.fasta  
- Plumapathes_pennacea_ORP_223_TMP.orthomerged.fasta  
- Hydractinia_polyclina_ORP_223_TMP.orthomerged.fasta  
- Stomolophus_meleagris_ORP_223_TMP.orthomerged.fasta  


Also in this directory are  
- the three "transcriptomes" of the cerianthid: hypostome_v2.ORP.fasta, body_v3.ORP.fasta, tentacle_v3.ORP.fasta (honestly I have no idea what the quality of these is, because almost all the reverse reads failed, so it is likely bad)  
- the Cerianthus genome as it now stands: ceri_pol5.fasta  
- the maker files: maker_opts.ctl (contains the controls for the maker run) and maker_evm.ctl, maker_bopts.ctl, and maker_exe.ctl (which just need to be in the same directory)  
- a text file containing a replica of the maker_opts.ctl file, in order to preserve the settings for this maker run for posterity: maker_opts1.txt  

Then I can run it using this code  
`module purge`  
`module load anaconda/colsa`  
`source activate maker-3.01.02`  

`cd /mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/`  

`export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/shared/augustus_config/config`  

`mpiexec -n 48 /mnt/lustre/macmaneslab/macmanes/test/maker/bin/maker -fix_nucleotides -base ceri_maker1 --ignore_nfs_tmp`  


During the run to check on the progress, I can use this code to see how many proteins it has so far  
`fasta_merge -d ceri_maker1.maker.output/ceri_maker1_master_datastore_index.log -o ceri_maker1`  

`grep -c ">" ceri_maker1.all.maker.proteins.fasta`  
Number of proteins as of 9-9-20: 17,124   
Number of proteins as of 9-14-20: 19,236  
Number of proteins as of 9-18-20: 22,841
Number of proteins at the end: 29,554  (Woo!)  

When this is finished (it took a long time because I put in basically every Cnid I could find), I ran the lines of code below in a script called "maker_end.sh". This is what generates all the output files I'll need to evaluate the run, etc. The first one is the one I used to check on the number of proteins as the run went along, but there are many more.   
`fasta_merge -d ceri_maker1.maker.output/ceri_maker1_master_datastore_index.log -o ceri_maker1`
`gff3_merge -d ceri_maker1.maker.output/ceri_maker1_master_datastore_index.log -o ceri_maker1.gff3 -n`
`sed -i '/^>/! s/-/N/g' ceri_maker1.all.maker.proteins.fasta`
`lastal -P38 /mnt/lustre/macmaneslab/macmanes/transporters/swissprot ceri_maker1.all.maker.proteins.fasta -f BlastTab > blast.out`
`maker_functional_fasta /mnt/lustre/macmaneslab/macmanes/transporters/uniprot_sprot.fasta blast.out` `ceri_maker1.all.maker.proteins.fasta > ceri_maker1.functional.proteins.fasta`
`maker_functional_fasta /mnt/lustre/macmaneslab/macmanes/transporters/uniprot_sprot.fasta blast.out ceri_maker1.all.maker.transcripts.fasta > ceri_maker1.functional.transcripts.fasta`
`maker_functional_gff /mnt/lustre/macmaneslab/macmanes/transporters/uniprot_sprot.fasta blast.out ceri_maker1.gff3 > ceri_maker1.functional.gff3`
`maker_map_ids --prefix Cerianthus_borealis --justify 6 ceri_maker1.gff3 > ceri_maker1.genome.all.id.map`
`map_fasta_ids ceri_maker1.genome.all.id.map  ceri_maker1.functional.proteins.fasta`
`map_gff_ids ceri_maker1.genome.all.id.map  ceri_maker1.functional.gff3`
`map_fasta_ids ceri_maker1.genome.all.id.map  ceri_maker1.functional.transcripts.fasta`


### Extras - for the second run of MAKER  

#### Long BUSCO run  

I am also going to do a run of busco with the "long" setting turned on, so that I can add this to the 2nd maker iteration. (And I finally got busco to run with multiple threads on a genome!)

This is the code I used:
`export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/config`

`run_BUSCO.py -i /mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/ceri_pol5.fasta -o ceri_pol5_long -m geno -l /mnt/lustre/hcgs/shared/databases/busco/metazoa_odb9 --long -c 40`  

I also have a file called "config.ini" in this same directory, and another directory called "config" that has multiple sub directories (that's where the augustus path leads). Finally got the right combination of things there, and off it goes (busco_long.sh).  

At the end, there is a directory that has a bunch of config files that maker will need. They are here: `/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/run_ceri_pol5_long/augustus_output/retraining_parameters/`  

I am making a new directory to hold these files in the maker_2020 directory I've been doing all of this in, and it's name will match the prefix of all the files in the augustus_output directory above. `/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/BUSCO_ceri_pol5_long_2574205824`. Then I can copy all of those files over to this new directory, and give maker the path to this new directory.  

#### Repeat Modeler    

I'll also include a run of RepeatModeler in the next run, code as follows:  
`BuildDatabase -name ceri ceri_pol5.fasta`  
`RepeatModeler -database ceri -pa 40`  

Results for this run can be found here: `/mnt/oldhome/macmaneslab/jlh1023/cerianthid/maker_2020/RM_172934.TueSep151359192020 ( consensi.fa.classified )`  


### Second run of MAKER  

Most things will stay the same for the second run of MAKER, but I will incorporate the two "extras" from above.  

Settings in the maker_opts.ctl file that I changed (these settings are preserved in maker_opts2.txt):  
- In the "Gene Prediction" section, put the path to the augustus files generated by the long busco run where it says "augustus_species": `augustus_species=/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/BUSCO_ceri_pol5_long_2574205824/`  
- In the "Repeat Masking" section, put the path to the repeat modeler output where it says "rmlib": `rmlib=/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/RM_172934.TueSep151359192020/consensi.fa.classified`  


### BUSCO scores for iterations of MAKER  

I'm only planning two MAKER runs, but if the busco scores are quite different between them, then I'll know it might make more sense to keep going.  

After I ran the ending stuff (right before "extras", above), I ran busco on the "ceri_maker1.all.maker.proteins.fasta" file, and it was unfortunately, a little lower than I'd hoped. The complete BUSCO score is 53.6% for the euk database, and 55.1% for the met one. Sadly, this is probably because the transcriptomic reads we have for this species are so bad and/or incomplete, and there isn't a ton we can do about that.  

### Adding in transcriptomes  

Paper just came out (Klompen, Macrander, Reitzel, and Stampar) that published four cerianthid transcriptomes, including ours! So I downloaded all the reads to a new directory (/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/transcriptomes).  

Then I ran the ORP on them using scripts like this "Isarachnanthus_nocturnus_orp.sh", and then I took the resulting transcriptomes, and popped them into the newest run of MAKER. The Pachycerianthus borealis one is in the "est" section, and the other three are in the "altest" section. I think having more closely related organisms and a better transcriptome is going to help a ton.  

I stopped the second run because I hadn't spread it across enough nodes, and it will have to run the first steps again too, which is unfortunate, but hopefully it won't take as long as last time.  

So the most up-to-date MAKER run is running from the script "maker_mpi.sh" in the same directory as all the others, but it will be spitting out files with "ceri_maker3" as the prefix, even though it's really more like 2.5.   
