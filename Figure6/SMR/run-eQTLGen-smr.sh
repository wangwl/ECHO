#!/bin/sh -l
# FILENAME: eQTLGen-smr.sh
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --mem=70g
#SBATCH --time=15:30:00

dir=/anvil/projects/x-mcb130189/Wubin/ECHO/rwang
smr=/home/x-rwang22/smr-1.3.1-linux-x86_64/smr

exposure=$1
output=$2
${smr} --beqtl-summary ${exposure} --beqtl-summary ${dir}/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --heidi-off --out ${output}

