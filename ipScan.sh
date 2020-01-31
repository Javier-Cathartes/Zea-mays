#!/usr/bin/env bash

#
#	ipScan
#	Performs a interproscan using PROSITE appl
#	:

file=$1
number_of_headers=$(grep ">" $file | wc -l)
number_of_lines=$(echo "2*$number_of_headers" | bc -l)
awk '{if ($1 ~ ">") {printf "\n"$0"\n"} else {printf $0}}' $file \
	| sed '/^$/d' > prot.tmp

output=${file}.tsv
> $output
for i in $(seq 40 40 $number_of_lines);
do
	print "                   ##"
	print "                   ##"
	print "BEGING PROC ipScan ##"
	print "#####################"
	print "# Processing sequences from $(echo "$i - 40 | bc -l) to $(echo "$i + 40" | bc -l) #"
	print "########################################"


	head -n $i prot.tmp | tail -n 40 > prot.launch
	interproscan -i prot.launch \
			-f tsv \
			-appl PROSITEPATTERNS,PROSITEPROFILES \
			--goterms
	cat prot.launch.tsv >> $output
	rm prot.launch.tsv
done
