# Bromeliads
Scripts for the work with data


# 1. Extraction of files 

bsub -e extraMinus.err -J ext_minus -o ext_minus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/minus_genes "tar -zxvf /data/ul/dee/dee_phylo/mserrano/Bromelia/FastaMinusGenes.tar.gz"

bsub -e extraPlus.err -J ext_mplus -o ext_plus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/plus_genes "tar -zxvf /data/ul/dee/dee_phylo/mserrano/Bromelia/FastaPlusGenes.tar.gz"


# 2. Counting N in alignments and size, then filtering 

##Estimating parameters 

bsub -e countingP.err -J Count_plus -o Count_plus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/6.fastaCorrectedPlusFinal/ "perl /scratch/beegfs/monthly/mserrano/all_genes/counting_nucleotides.pl > /scratch/beegfs/monthly/mserrano/all_genes/align_length_plus.out"

bsub -e countingM.err -J Count_minus -o Count_minus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/minus_genes/7.fastaCorrectedMinusRevCompFinal "perl /scratch/beegfs/monthly/mserrano/all_genes/counting_nucleotides.pl > /scratch/beegfs/monthly/mserrano/all_genes/align_length_minus.out"

# In R 

counts_plus <- read.table("align_length_plus.out", h=F, sep="\t")
counts_minus <- read.table("align_length_minus.out", h=F, sep="\t")

dim(counts_plus)
length( which(counts_plus$V2 > 300 & counts_plus$V3 < 0.5) )
length( which(counts_plus$V2 > 300 & counts_plus$V3 < 0.3) )

dim(counts_minus)
length( which(counts_minus$V2 > 300 & counts_minus$V3 < 0.5) )
length( which(counts_minus$V2 > 300 & counts_minus$V3 < 0.3) )

# Filtering files for the first length 

bsub -e filter_plus.err -J filter_plus -o filter_plus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/plus_genes/6.fastaCorrectedPlusFinal/ "bash /scratch/beegfs/monthly/mserrano/all_genes/first_filtering_files.sh"

bsub -e filter_minus.err -J filter_minus -o filter_minus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/minus_genes/7.fastaCorrectedMinusRevCompFinal/ "bash /scratch/beegfs/monthly/mserrano/all_genes/first_filtering_files.sh"


## This first filter stored the fasta files in : 
## sftp://mserrano@prd.vital-it.ch/data/ul/dee/dee_phylo/mserrano/Bromelia/all_genes/plus_genes/6.fastaCorrectedPlusFinal/filtered_alignments
## sftp://mserrano@prd.vital-it.ch/data/ul/dee/dee_phylo/mserrano/Bromelia/all_genes/minus_genes/7.fastaCorrectedMinusRevCompFinal/filtered_alignments

# Counting files 

ls -1 6.fastaCorrectedPlusFinal/filtered_alignments/ | wc -l
ls -1 7.fastaCorrectedMinusRevCompFinal/filtered_alignments/ | wc -l


###############################################################################
3. Looking for proteins with many stop codons, many Ns in the sequence,to remove them and filter files
###############################################################################

# converting all files into proteins to identify those with stop codons...before running the test of guidance...

bsub -e transl.err -J transl -o transl.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/test_alignments/ "perl translating_to_protein.pl"
bsub -e transl.err -J transl -o transl.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "perl translating_to_protein.pl"


# These file contains all the alignments converted into proteins for the codon-stop evaluation (but were erased by vital-it)

bsub -e zip_prots.err -J zip_prots -o zip_prots.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "tar -cjf transl_protein_plus_minus.tbz 6.fastaCorrectedPlus_protein/ 7.fastaCorrectedMinus_protein/"

bsub -e zip_nt_filtered_used.err -J zip_nt_filtered_used -o zip_nt_filtered_used.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "tar -cjf zip_nt_filtered_used.tbz 6.fastaCorrectedPlusFinal/filtered_alignments/ 7.fastaCorrectedMinusRevCompFinal/filtered_alignments/"

## checking for stop codons in all the protein files 

grep -rl 7.fastaCorrectedMinus_protein/ -e "*" > stop_codons_proteins_minus.txt 
grep -rc 7.fastaCorrectedMinus_protein/ -e "*" > counts_stop_codons_proteins_minus.txt

wc -l stop_codons_proteins_minus.txt 
grep -rl 6.fastaCorrectedPlus_protein/ -e "*" > stop_codons_proteins_plus.txt
wc -l stop_codons_proteins_plus.txt 

## checking in R

counts_plus <- read.table("counts_stop_codons_proteins_plus.txt",sep=":")
genes_plus_with_stop <- counts_plus [ which(counts_plus[,2] > 0) , ]
length(which(genes_plus_with_stop[,2]==1))

counts_minus <- read.table("counts_stop_codons_proteins_minus.txt",sep=":")
genes_minus_with_stop <- counts_minus [ which(counts_minus[,2] > 0) , ]
length(which(genes_minus_with_stop[,2]==1))

pdf("Number_stopCodons_genes.pdf", paper="a4r")
par(mfrow=c(2,1))
barplot(table(genes_plus_with_stop[,2]), xlab="# stop codons", ylab="# genes")
barplot(table(genes_minus_with_stop[,2]), xlab="# stop codons", ylab="# genes")
dev.off()

length(which(genes_plus_with_stop[,2]>29))
length(which(genes_minus_with_stop[,2]>29))


# with the script checking_evaluation.R we selected the genes: 

 where maximum 2 seqs have > 90% bad sequence ==N then they were logged in file "genes_to_rem_max2seq_minus.txt"
 where 3 or more seqs have >90% bad sequence == N then all the gene was removed, list genes in file "genes_to_rem_allGene_minus.txt"
 
 # these files are for plus and minus 
 # files genes_to_rem_allGene_* used to run the checking_codon_length_divided3_%Ns.R and discard those genes 
 
## With the this file generate which alignments pass the filter 

checking_codon_length_divided3_%Ns.R 

producing this files: 	files_to_keep_plus_d3_stop29_filter.txt
						files_to_keep_minus_d3_stop29_filter.txt

## With all these criteria we proceed to next filter which is the non-divided by 3 and more than 29 stop codons 
## sending the filtering in all files 

bsub -e filter134_plus.err -J filter134_plus -o filter134_plus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "bash /scratch/beegfs/monthly/mserrano/all_genes/filtering_files_all_filters.sh"
 
bsub -e filter134_minus.err -J filter134_minus -o filter134_minus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "bash /scratch/beegfs/monthly/mserrano/all_genes/filtering_files_all_filters_minus.sh"

 
 ## now the individual sequences have to be removed from the alignment 
 
 # use the file removing_Nseqs.sh and the file genes_to_rem_max2seq_plus.txt
 # sending the script for removing the one or two sequences that are >90% N 

bsub -e rem_seqs_plus.err -J rem_seqs_plus -o rem_seqs_plus.out -q normal -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "bash /scratch/beegfs/monthly/mserrano/all_genes/removing_Nseqs.sh"

bsub -e rem_seqs_minus.err -J rem_seqs_minus -o rem_seqs_minus.out -q normal -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "bash /scratch/beegfs/monthly/mserrano/all_genes/removing_Nseqs.sh"


#####################################
3. Creating the phylogenetic trees 
#####################################

# For all the files that were not modified in plus and minus we execute phyml as follows (scripts 9-13): (/scratch/beegfs/monthly/mserrano/all_genes/plus_genes/6.fastaCorrectedPlusFinal/filtered_alignments)

# Warning always remember to change the path to the MSA alignments in both files 
# use only this version otherwise the trees complain 
module add Phylogeny/phyml/3.3.20170119;

# local 
# /home/marthaserrano/Downloads/PhyML-3.1/PhyML-3.1_linux32 -i Tillandsia_29samples_UG_Aco001189_minDP3.phy -q -m HKY85 -b -1 -o tlr -o tlr --constraint_file RAxML_bipartitions_constraint.tre -u RAxML_bipartitions_constraint.tre

## check that constrained tree exists in the folder it should RAxML_bipartitions_constraint.tre

python 04_writePhyml_LSFscript_andMap_BROMELIA.py 04.1_phyml_childprocess_BROMELIA.py /scratch/beegfs/monthly/mserrano/all_genes/minus_genes/7.fastaCorrectedMinusRevCompFinal/ dee-long
 

python 04_writePhyml_LSFscript_andMap_BROMELIA.py 04.1_phyml_childprocess_BROMELIA.py /scratch/beegfs/monthly/dscherre/plus_genes/6.fastaCorrectedPlusFinal/ normal

bsub < RUN_phyml_lsf_array.sh


## for the alignments that were modified (in new_phyml_runs/ folder) it requires that the constrain tree is re-written because not all taxa are included 


export PATH=/software/bin:$PATH;
module add R/3.0.2;

bsub -e const_plus.err -J const_plus -o const_plus.out -q dee-long -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "R --no-save < ./parsing_align_outgroup_constraints.R"

bsub -e const_minus.err -J const_minus -o const_minus.out -q dee-hugemem -cwd /scratch/beegfs/monthly/mserrano/all_genes/ "R --no-save < ./parsing_align_outgroup_constraints.R"

sftp://mserrano@prd.vital-it.ch/scratch/beegfs/monthly/mserrano/all_genes/new_phyml_runs

python 04_writePhyml_LSFscript_andMap_BROMELIA_modifruns.py 04.1_phyml_childprocess_BROMELIA_modifruns.py /scratch/beegfs/monthly/mserrano/all_genes/new_phyml_runs/ dee-hugemem

bsub < RUN_phyml_lsf_array.sh
