#!/bin/sh

cd /scratch1/zheyuli/indiv_gt

split -l 200000 -a 3 -d --additional-suffix=.vcf body.vcf --verbose

