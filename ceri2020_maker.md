## Running MAKER on the new Cerianthus assembly

I need to run MAKER on my new assembly, so first I am gathering all of the other Cnidarians that I will use in the run. All of these will be in this directory: /mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/  

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
`module purge  
module load anaconda/colsa  
source activate maker-3.01.02`  

`cd /mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/`  

`export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/shared/augustus_config/config`  

`mpiexec -n 48 /mnt/lustre/macmaneslab/macmanes/test/maker/bin/maker \  
-fix_nucleotides -base ceri_maker1 --ignore_nfs_tmp`  


During the run to check on the progress, I can use this code to see how many proteins it has so far  
`fasta_merge -d ceri_maker1.maker.output/ceri_maker1_master_datastore_index.log -o ceri_maker1`  

`grep -c ">" ceri_maker1.all.maker.proteins.fasta`
Number of proteins as of 9-9-20: 17,124  
Number of proteins as of 9-14-20: 19,236  


### Extras  

#### Long BUSCO run  

I am also going to do a run of busco with the "long" setting turned on, so that I can add this to the 2nd maker iteration. (And I finally got busco to run with multiple threads on a genome!)

This is the code I used:
`export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/config`

`run_BUSCO.py -i /mnt/lustre/macmaneslab/jlh1023/cerianthid/maker_2020/ceri_pol5.fasta \
-o ceri_pol5_long -m geno -l /mnt/lustre/hcgs/shared/databases/busco/metazoa_odb9 --long -c 40`  

I also have a file called "config.ini" in this same directory, and another directory called "config" that has multiple sub directories (that's where the augustus path leads). Finally got the right combination of things there, and off it goes (busco_long.sh).  

#### Repeat Modeler    

I'll also include a run of RepeatModeler in the next run, code as follows:  
`BuildDatabase -name ceri ceri_pol7.fasta  
RepeatModeler -database ceri -pa 40`  

Results for this run can be found here: `/mnt/oldhome/macmaneslab/jlh1023/cerianthid/maker_2020/RM_172934.TueSep151359192020 ( consensi.fa.classified )`
