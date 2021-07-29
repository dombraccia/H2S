#!/bin/bash

## a script for grabing up to the family level of all hits from GutFunFind
## and pairing them with their gene hit info from {function}.genes.tsv
## (for supplementary table in H2S paper)

rm -f results/from-GutFunFind/Cysteine_Degradation.genes_table.txt
rm -f results/from-GutFunFind/family_level_lineage_tmp.txt
rm -f results/from-GutFunFind/genes_modified_tmp.txt

cut -d ';' -f5-7 results/from-GutFunFind/Cysteine_Degradation.taxa_hits.txt >> results/from-GutFunFind/family_level_lineage_tmp.txt

genes=$(grep -v 'absent' results/from-GutFunFind/Cysteine_Degradation.genes.tsv | cut -f2)
for line in $genes; do
	#echo "starting line: $line"
	line=$(echo $line | sed -e $'s/,/\\\n/g' | sort -u | tr '\n' ',' | sed 's/.$//')
	#echo "ending line: $line"
	printf "%b" "$line\n" >> results/from-GutFunFind/genes_modified_tmp.txt
done

## combining lineage and processed gene columns
#genes=$(cat results/from-GutFunFind/genes_modified_tmp.txt)
paste results/from-GutFunFind/family_level_lineage_tmp.txt results/from-GutFunFind/genes_modified_tmp.txt > tmp.txt
sed '1 s/^/Bacterial Species\tGene hits to species\n/' tmp.txt > results/from-GutFunFind/Cysteine_Degradation.genes_table.txt

rm -f results/from-GutFunFind/genes_modified_tmp.txt
rm -f results/from-GutFunFind/family_level_lineage_tmp.txt
rm -f tmp.txt 

