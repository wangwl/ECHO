# Preprocessing
coloc_preprocessing.ipynb

The processed gwas catalog data were uploaded onto figshare ([https://figshare.com/articles/dataset/gwas_catalog_202406/26210795](https://figshare.com/articles/dataset/gwas_catalog_202406/26210795))

# Run coloc
```shell
outdir="/home/x-wding2/project_wubin/ECHO/coloc_results"
mkdir -p ${outdir}
rm -f run_coloc.sh

# cc.no_beta
for file in `ls ~/project_wubin/ECHO/cis-meQTL-norminal/*`; do
    bname=$(basename ${file})
    name=${bname/.tsv/}
    echo ${name}
    if ! [ -f ${outdir}/${name}.cc.no_beta.ipynb ]; then
        echo "papermill -p infile1 ${file} -p infile2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_cc_no_beta.tsv.gz -p phe1 dmr_id -p phe2 gwas_id -p has_beta1 no -p has_beta2 no -p N1 110 -p meta2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_harmonised_metadata.tsv -p outdir ${outdir} -p name ${name}.cc.no_beta -k r4.3.3 run_coloc_quant_cc.ipynb ${outdir}/${name}.cc.no_beta.ipynb" >> run_coloc.sh
    fi
done;

# cc.with_beta
for file in `ls ~/project_wubin/ECHO/cis-meQTL-norminal/*`; do
    bname=$(basename ${file})
    name=${bname/.tsv/}
    echo ${name}
    if ! [ -f ${outdir}/${name}.cc.with_beta.ipynb ]; then
        echo "papermill -p infile1 ${file} -p infile2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_cc_with_beta.tsv.gz -p phe1 dmr_id -p phe2 gwas_id -p has_beta1 no -p has_beta2 yes -p N1 110 -p meta2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_harmonised_metadata.tsv -p outdir ${outdir} -p name ${name}.cc.with_beta -k r4.3.3 run_coloc_quant_cc.ipynb ${outdir}/${name}.cc.with_beta.ipynb" >> run_coloc.sh
    fi
done;

# cc.no_beta
for file in `ls ~/project_wubin/ECHO/cis-meQTL-norminal/*`; do
    bname=$(basename ${file})
    name=${bname/.tsv/}
    echo ${name}
    if ! [ -f ${outdir}/${name}.quant.no_beta.ipynb ]; then
        echo "papermill -p infile1 ${file} -p infile2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_quant_no_beta.tsv.gz -p phe1 dmr_id -p phe2 gwas_id -p has_beta1 no -p has_beta2 no -p N1 110 -p meta2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_harmonised_metadata.tsv -p outdir ${outdir} -p name ${name}.quant.no_beta -k r4.3.3 run_coloc_quant_quant.ipynb ${outdir}/${name}.quant.no_beta.ipynb" >> run_coloc.sh
    fi
done;

# cc.with_beta
for file in `ls ~/project_wubin/ECHO/cis-meQTL-norminal/*`; do
    bname=$(basename ${file})
    name=${bname/.tsv/}
    echo ${name}
    if ! [ -f ${outdir}/${name}.quant.with_beta.ipynb ]; then
        echo "papermill -p infile1 ${file} -p infile2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_quant_with_beta.tsv.gz -p phe1 dmr_id -p phe2 gwas_id -p has_beta1 no -p has_beta2 yes -p N1 110 -p meta2 /home/x-wding2/project_wubin/ECHO/gwas_catalog/gwas_catalog_harmonised_metadata.tsv -p outdir ${outdir} -p name ${name}.quant.with_beta -k r4.3.3 run_coloc_quant_quant.ipynb ${outdir}/${name}.quant.with_beta.ipynb" >> run_coloc.sh
    fi
done;
```

## merge coloc results

```python
import os.path
import glob
import numpy as np

indir = os.path.expanduser("~/project_wubin/ECHO/coloc_results")
R=[]
for infile in glob.glob(os.path.join(indir, "*.txt")):
    bname=os.path.basename(infile)
    name = bname.split('.')[0]
    df = pd.read_csv(infile, sep='\t')
    df['Name']=name
    df['file']=bname
    R.append(df)
df=pd.concat(R)
df.sort_values('PP4',ascending=False,inplace=True)
# annotate gwas
df_gwas=pd.read_csv(os.path.expanduser("~/project_wubin/ECHO/gwas_catalog/gwas_catalog_harmonised_metadata.tsv"),
                    sep='\t',index_col='gwas_id')
for col in ['trait_description','ontology_mapping','sample_size','n_case','n_control']:
    df[col]=df.phe2.map(df_gwas[col].to_dict())
df.trait_description=df.trait_description.apply(lambda x:';'.join(eval(x)) if not pd.isna(x) else np.nan)
df.ontology_mapping=df.ontology_mapping.apply(lambda x:';'.join(eval(x)) if not pd.isna(x) else np.nan)
df.to_csv("coloc_results.tsv",sep='\t')
df.loc[(df.PP4>=0.8) & (df.n_snps >= 20)]
```

# Visualization
coloc_visualization.ipynb

plot.ipynb