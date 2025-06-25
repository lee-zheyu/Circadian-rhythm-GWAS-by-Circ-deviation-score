#!/bin/sh

cd ~/cmb0/circadian_gene/data/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/indiv_gt

tail -n +84 indiv_gt.vcf > body.vcf

