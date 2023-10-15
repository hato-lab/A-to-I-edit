#!/usr/bin/env bash

## to check A-to-I edit at pre-specified sites
## this script was used for sanity check

### this configuration is sequence kit specific (F2R1/F1R2 and FLAG) 
is_positive_strand_gene="positive" 
geneName="Khdc4"
strandedness="F2R1"
flag_1="163"
flag_2="83"
my_chr="chr3"
my_position="88704183"


my_directory="~/AtoI"
my_reference="GRCm38.primary_assembly.genome.fa"

mkdir ${my_directory}/${geneName}
cd ${my_directory}/BAM

awk -v is_positive_strand_gene=${is_positive_strand_gene} -v geneName=${geneName} -v strandedness=${strandedness} -v flag_1=${flag_1} -v flag_2=${flag_2} -v my_chr=${my_chr} -v my_position=${my_position} \
'BEGIN{FS="\t";OFS="\t"} FNR == 1{print is_positive_strand_gene; print geneName; print strandedness; print flag_1; print flag_2; print my_chr; print my_position}' `echo ${my_directory}/dummy.txt` > ${my_directory}/search_condition.txt


array1=(`ls | egrep bam$`)

### make -f flag specified bam files, concatenate, sort, and index 

for i in ${array1[@]}; do \
	samtools view -h -f ${flag_1} ${i} ${my_chr}:${my_position}-${my_position} | awk 'BEGIN{FS="\t";OFS="\t"} {print $0}' >> ${geneName}_${i}; \
	samtools view -f ${flag_2} ${i} ${my_chr}:${my_position}-${my_position} | awk 'BEGIN{FS="\t";OFS="\t"} {print $0}' >> ${geneName}_${i}; \
	samtools view -S -b ${geneName}_${i} > ${strandedness}_${geneName}_${i}; \
	samtools sort -o ${strandedness}_${geneName}_${i}_sorted.bam ${strandedness}_${geneName}_${i} ; \
	samtools index ${strandedness}_${geneName}_${i}_sorted.bam ; \
	awk -v filename=${i} 'BEGIN{FS="\t";OFS="\t"} FNR == 1{print "@", filename}' `echo ${strandedness}_${geneName}_${i}_sorted.bam` >> foo.tsv ; \
	bam-readcount -w1 -f ${my_reference} ${strandedness}_${geneName}_${i}_sorted.bam ${my_chr}:${my_position}-${my_position} >> foo.tsv ; \
	mv ${strandedness}_${geneName}_${i}_sorted.bam ${my_directory}/${geneName}/ ; \
	mv ${strandedness}_${geneName}_${i}_sorted.bam.bai ${my_directory}/${geneName}/ ; \
	rm ${geneName}_${i} ; \
	rm ${strandedness}_${geneName}_${i}; done


mv foo.tsv ${my_directory}/${geneName}/


cd ${my_directory}/${geneName}
awk -v strand=${strandedness} -v chr=${my_chr} -v position=${my_position} 'BEGIN{FS="\t";OFS="\t"} \
FNR == NR && $1 ~ /^@/{filename_array[FNR]=$2; next} \
FNR != NR && $1 ~ chr && $2 ~ position && length($2) == length(position) { \
	A_split=split($6, A_array, ":"); C_split=split($7, C_array, ":"); G_split=split($8, G_array, ":"); T_split=split($9, T_array, ":"); \
	if (strand == "F1R2" && (C_array[2] + T_array[2]) > 0){ AtoI = (C_array[2]/(C_array[2] + T_array[2])) } \
	if (strand == "F1R2" && (C_array[2] + T_array[2]) == 0){ AtoI = 0 } \
	if (strand == "F2R1" && (A_array[2] + G_array[2]) > 0){ AtoI = (G_array[2]/(A_array[2] + G_array[2])) }; \
	if (strand == "F2R1" && (A_array[2] + G_array[2]) == 0){ AtoI = 0 }; \
	print filename_array[(FNR-1)], $1, $2, "ref", $3, $4, A_array[1], A_array[2], C_array[1], C_array[2], G_array[1], G_array[2], T_array[1], T_array[2], strand, "AtoI%", AtoI}' \
	foo.tsv foo.tsv > foo_final_${geneName}_${my_chr}_${my_position}_${strandedness}.tsv
	
