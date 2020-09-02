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
