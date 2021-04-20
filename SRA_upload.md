## Uploading the cerianthid genome to the SRA  

Since Dave has already submitted the Illumina reads to the SRA, and the nanopore fast5 files are being found, I am going to focus on getting the genome and its annotation into the correct format. Happening in this directory: /mnt/lustre/macmaneslab/jlh1023/cerianthid/submit_sra/  

There are about 40 different sites that NCBI has put together to tell you how to do this, and each only has partial information that you need. Some of these are below. Stay strong.  
- https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/  
- https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/  
- https://submit.ncbi.nlm.nih.gov/genbank/template/submission/  
- https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/  
- https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/  

This last site has the template that I used to make the command I used to do the conversion. But first you have to make sure that table2asn is installed (note: this is a different thing than tbl2asn, so do not get them confused.). It can be installed from here: https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/, but I'll be honest, I got Toni to do it for me because this submission was already a **whole thing**.  

Another thing to note is that even though the guides start off with the command `table2asn`, do not be fooled! The command is actually `table2asn_GFF`, just to keep things confusing. The command I used is below.  
`table2asn_GFF -M n -J -c w -euk -t template.sbt -gaps-min 10 -l paired-ends -j "[organism=Pachycerianthus borealis]" -i Pachycerianthus_borealis.fa -f pachyceri_maker.functional.gff3 -o Pachycerianthus_borealis.sqn -Z -locus-tag-prefix JQW98`  

It looks like this worked, so I sent the path to the new file (/mnt/lustre/macmaneslab/jlh1023/cerianthid/submit_sra/Pachycerianthus_borealis.sqn, although it did generate other files too) to Dave, and he should be able to upload it to the bioproject with the rest of the info from the same individual.  


**Update**  
The genome assembly didn't get accepted right away because there is a sequence in it that should not be there. So I'm going to remove that sequence, make sure the rest of the genome sequences meet the requirements of genbank, and remake the file.  

The offending sequence is Pachycerianthus_borealis_2734  
I unwrapped the fasta first, since it previously had annoying linefeeds in all the sequences.  
`fasta_formatter -w 0 -i Pachycerianthus_borealis.fa -o Pachycerianthus_borealis_unwrap.fa`  

Then I removed the sequence.  
`sed "/Pachycerianthus_borealis_2734/,+1d" Pachycerianthus_borealis_unwrap.fa > Pachycerianthus_borealis_clean.fa`  

To confirm:  
`wc -l Pachycerianthus_borealis_unwrap.fa`  11666  
`wc -l Pachycerianthus_borealis_clean.fa`  11664  

Looks like it got removed, and when I grep for the sequence name it no longer returns any results.  

They also want me to make sure that there aren't a ton of Ns at the ends of sequences, and that they aren't too short. I know we don't have any short seqs in here, and I grepped for "N" and didn't find any matches, so I think we are good there. Now I'll just rerun the command above to remake the sqn file. I renamed the old sqn file to "Pachycerianthus_borealis_v1.sqn" so that I can keep the same output file name.     
`table2asn_GFF -M n -J -c w -euk -t template.sbt -gaps-min 10 -l paired-ends -j "[organism=Pachycerianthus borealis]" -i Pachycerianthus_borealis_clean.fa -f pachyceri_maker.functional.gff3 -o Pachycerianthus_borealis.sqn -Z -locus-tag-prefix JQW98`  

Note that this process is submitted as a job with this script: /mnt/oldhome/macmaneslab/jlh1023/cerianthid/submit_sra/conversion.sh 

