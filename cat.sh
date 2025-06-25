#!/bin/sh

cd /scratch1/zheyuli/indiv_gt

for ((i=0; i<=9; i++)); do
	inputfile="x00$i.vcf"
	outputfile="p00$i.vcf"
	cat head.vcf $inputfile > $outputfile
done


for ((i=10; i<=99; i++)); do
	inputfile="x0$i.vcf"
	outputfile="p0$i.vcf"
	cat head.vcf $inputfile > $outputfile
done

for ((i=100; i<=348; i++)); do
	inputfile="x$i.vcf"
	outputfile="p$i.vcf"
	cat head.vcf $inputfile > $outputfile
done



