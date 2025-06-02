
module load r/4.0.5

scripts=/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/smr/scripts
eQTLGen=/anvil/projects/x-mcb130189/Wubin/ECHO/rwang/eQTLGen
smr=/home/x-rwang22/smr-1.3.1-linux-x86_64/smr

sbatch ${scripts}/call-find-snps.sh ${eQTLGen}/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi
${smr} --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --update-esi cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi.snp.cor.txt

sbatch ${scripts}/call-find-genes.sh ${eQTLGen}/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi
${smr} --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --update-epi cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi.gene.cor.txt

# dir=/anvil/projects/x-mcb130189/Wubin/ECHO/rwang
# ${smr} --beqtl-summary ${dir}/smr/cis-meQTL/B-Naive.permutation --beqtl-summary ${dir}/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --heidi-off --out nB_eQTLGen_
