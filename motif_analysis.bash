#!/usr/bin/env bash


## extract sequences around A-to-I edit sites

masterfile="masterfile_above_20_counts_in_more_than_3_samples.tsv"
filter_name="hyperedit_pos_filter_above20"

awk 'BEGIN{FS="\t"; OFS="\t"} FNR == 1{print $0} FNR != 1 && $82 ~ /TRUE/{print $0}' ${masterfile} > ${masterfile}_${filter_name}.tsv


Edit_list=${masterfile}_${filter_name}.tsv
subset1="only_3UTR"

mkdir ${filter_name}_and_${subset1}_motif
cd ${filter_name}_and_${subset1}_motif


###	$39 = genomicLocation 							
### extract 50 nt +/- from the edit site                         				
### split files into + and - reads	(redirecting inside of awk command)	
### start position -1 (UCSC is 0-based)

awk -v subset1=${subset1} 'BEGIN{FS="\t"; OFS="\t"} FNR == 1{next} \
	FNR != 1 && subset1 != /all/ && $39 ~ subset1{split($1, new_array, " "); \
		if (new_array[3] ~ /+/){ \
			edit_upstream = (new_array[2] - 50); edit_downstream = (new_array[2] + 50); \
			print new_array[1] ":" (edit_upstream - 1) "-" edit_downstream > "for_twoBitToFa_pos_strand_input.txt" }; \
		if (new_array[3] ~ /-/){ \
			edit_upstream = (new_array[2] + 50); edit_downstream = (new_array[2] - 50); \
			print new_array[1] ":" (edit_downstream -1) "-" edit_upstream > "for_twoBitToFa_neg_strand_input.txt" }}' ${Edit_list}



~/bin/twoBitToFa mm10.2bit ${subset1}_pos_strand_genes.fa -seqList="for_twoBitToFa_pos_strand_input.txt" -noMask

~/bin/twoBitToFa mm10.2bit ${subset1}_neg_strand_genes.fa -seqList="for_twoBitToFa_neg_strand_input.txt" -noMask

## reverse complement the neg strand
~/bin/faRc ${subset1}_neg_strand_genes.fa ${subset1}_neg_strand_genes_RC.fa -keepName

## flatten the fasta. Also conver to upper for neg strand genes
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $0; next} $0 !~ /^>/{printf "%s", $0} $0 ~ /^>/{printf "\n%s\n", $0} END{printf "\n"}' ${subset1}_pos_strand_genes.fa > ${subset1}_pos_strand_genes_flat.fa
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{print $0; next} $0 !~ /^>/{printf "%s", toupper($0)} $0 ~ /^>/{printf "\n%s\n", $0} END{printf "\n"}' ${subset1}_neg_strand_genes_RC.fa > ${subset1}_neg_strand_genes_RC_flat.fa


### combine pos and neg fa files 
cat ${subset1}_pos_strand_genes_flat.fa ${subset1}_neg_strand_genes_RC_flat.fa > ${subset1}_final.fa

~/bin/kpLogo ${subset1}_final.fa

