# Cell Ranger
* `Cell Ranger` is a tool for converting fastq outputs of a sequencer to expression matrix file avilable for further analysis in R and Python.
```bash
cellranger count \
	--id=HRR864929 \
	--transcriptome=/home/flynn/ref/refdata-gex-GRCh38-2020-A \
	--fastqs=/media/flynn/RuochongSSD/gsa/HRR864929/R_0 \
	--sample=R21058446-OES211029150-OES211029150 \
	--create-bam=true \ # bam file for further analysis
	--nosecondary  \ # defualt analysis provided by Cell Ranger is not needed
	--output-dir /home/flynn/HRR864929 \
	--nopreflight \ # no need to check Cell Ranger working environment 
  	--disable-ui # no need to monitor Cell Ranger performance
```

# Velocyto.py
* `Velocyto` is necessary to perform RNA Velocity Trajectory Analysis. Further data analysis and visulization is available in Python.
```bash
conda activate velocyto
velocyto run \
	-b "/home/flynn/HRR864929/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" \
	-o "/home/flynn/HRR864929/outs/loom"\
	 -m "/home/flynn/ref/repeat_msk.gtf" \
	 "/home/flynn/HRR864929/outs/possorted_genome_bam.bam" \
	 "/home/flynn/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
```
