#!/bin/bash

########################################################
# Script example to simplify the use of many options #
########################################################

#PATH to the FASTA file used to produce the annotation
GENOME=`dirname "$0"`"/Pachycerianthus_borealis.fa"

#PATH to the ANNOTATION in gff3 FORMAT
ANNOTATION=`dirname "$0"`"/pachyceri_maker.functional.gff3"

#PROJECT name registered on EMBL
PROJECT="PRJEB42842 "

#Locus tag registered on EMBL
LOCUS_TAG="PBOREALIS"

# species name
SPECIES="2736680"

# Taxonomy
#TAXONOMY="INV"

#The working groups/consortia that produced the record. No default value
#REFERENCE_GROUP="XXX"

#Translation table
TABLE="1"

#Molecule type of the sample.
MOLECULE="genomic DNA"

myCommand="EMBLmyGFF3 -i $LOCUS_TAG -p $PROJECT -m \"$MOLECULE\" -r $TABLE -t linear -s \"$SPECIES\" -o Pachycerianthus_borealis_genomeassembly_1.0.embl $ANNOTATION $GENOME $@"
echo -e "Running the following command:\n$myCommand"

#execute the command
eval $myCommand
