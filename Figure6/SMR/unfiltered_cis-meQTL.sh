smr=/home/x-rwang22/smr-1.3.1-linux-x86_64/smr
type=$1

###### myeqtl.esi file - SNPs ######
# A1 = Ref 
# A2 = Alt
# awk 'BEGIN{OFS="\t"} NR!=1 {sub(/chr/, "",$1); print $1, $3, 0, $2, $6, $7}' ${merged} > cis-meQTL.esi

cd /anvil/projects/x-mcb130189/Wubin/ECHO/rwang/smr/cis-meQTL

[ -d ${type} ] || mkdir -p "${type}"

cd ${type}

###### myeqtl.epi file - probes ######
# dmrs=/anvil/projects/x-mcb130189/Wubin/ECHO/permute_cis.gDMR.bed 
# cd /anvil/projects/x-mcb130189/Wubin/ECHO/rwang/smr/cis-meQTL

# create a probe file for each cell type B-Naive.tsv 
# awk -v var="${type}" 'BEGIN{OFS="\t"} $5 == var {sub(/chr/, "",$1); print $1, $4, 0, $2, $4, "+"}' ${dmrs} | sort | uniq > ${type}_dmr.epi

# make the besd file
meQTL=/anvil/projects/x-mcb130189/Wubin/ECHO/raw.cis-meQTL
input=/anvil/scratch/x-wangwl/02.mQTL/04.bulkQTL/04.cis/02.permute

# process the besd files first: 
zcat ${input}/${type}.permutation.gz | awk 'BEGIN{OFS="\t"} NR != 1 && $8 != "NA" {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $18, $19, $20}' > ${type}.permutation.filtered.smr.tsv

# process the besd files first: 
zcat ${input}/${type}.permutation.gz | awk 'BEGIN{OFS="\t"} $8 != "NA" {sub(/chr/, "",$2); print $2, $1, 0, $3, $1, "+"}' | sort | uniq > ${type}_dmr.epi

${smr} --eqtl-summary ${type}.permutation.filtered.smr.tsv --qtltools-permu-format --make-besd --out ${type}.permutation
${smr} --beqtl-summary ${type}.permutation --update-epi ${type}_dmr.epi
${smr} --beqtl-summary ${type}.permutation --update-esi ../cis-meQTL.esi


