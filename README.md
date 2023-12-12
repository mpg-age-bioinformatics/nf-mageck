# nf-fastqc

Create the test directory:
```
mkdir -p ~/nf-kallisto-test/raw_data
```

Download the demo data:
```
cd ~/nf-kallisto-test/raw_data
curl -J -O https://datashare.mpcdf.mpg.de/s/jcEaS5vqpJO0lOy/download
curl -J -O https://datashare.mpcdf.mpg.de/s/XHanbnjfvQ9rACD/download
curl -J -O https://datashare.mpcdf.mpg.de/s/sIebkRdMfMSweq2/download
curl -J -O https://datashare.mpcdf.mpg.de/s/zoNxS9vRI7jl77y/download
curl -J -O https://datashare.mpcdf.mpg.de/s/0WHGNIhjJC792lY/download
curl -J -O https://datashare.mpcdf.mpg.de/s/ZlM0lWKPh8KrP6B/download
curl -J -O https://datashare.mpcdf.mpg.de/s/o3O6BKaEXqB7TTo/download
```

Download the paramaters file:
```
cd ~/nf-kallisto-test
curl -J -O https://raw.githubusercontent.com/mpg-age-bioinformatics/nf-fastqc/main/params.json
```

Run the workflow:
```
RELEASE=1.0.0
ORIGIN="mpg-age-bioinformatics/"

nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry images && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry pre_process && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_count && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_pretest && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_test && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_pathway && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_plot && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_premle && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_mle && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_vispr && \
nextflow run ${ORIGIN}nf-mageck ${RELEASE} -params-file params.json -entry mageck_flute
```

## Parameters

```
{ 
  # where should the results be stored
  "project_folder" : "/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data" ,
  
  # an excel file with the library sheet, the sample and sampleNames sheets
  "reference_file":"/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/20231012.092020.f_xfkxly.crispr.xlsx",
  
  # raw data folder
  "raw_fastq":"/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/Omizt65/",
  
  # where the samples with meaningfull names should be stored after relabeling (ie. ln -s during preprocessing)
  "renamed_fastq":"/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/raw_renamed",
  
  # fastqc input
  "fastqc_raw_data" : "/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/raw_renamed" ,
  
  # cutadapt input
  "cutadapt_raw_data" : "/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/raw_renamed",
  
  # sgRNA_size
  "sgRNA_size" : "30",
  
  # sgRNA size used for calculating sgRNA efficiencies
  "SSC_sgRNA_size" : "20",
  
  # upstream sequence
  "upstreamseq" : "TGTGGAAAGGACGAAACACC",
  
  # path to library.csv as generated during preprocessing
  "library" : "/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/library.csv",
  
  # mageck count input 
  "input_count" : "/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/cutadapt_output" , 
  
  # mageck count output in project_folder
  "output_count" : "mageck_output/fastq",
  
  # how the counts file was generated, default: mageck ie. from fastq files processed in mageck count
  "mageck_counts_type":"mageck",
  
  # samples.tsv as generated during preprocessing
  "samples_tsv":"/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/samples.tsv",
  
  # mageck count output in project_folder in project_folder
  "output_test":"mageck_output/fastq/test",
  
  # Whether to remove zero-count sgRNAs in control and/or treatment experiments (none, control, treatment, both, any)
  "mageck_test_remove_zero": "none",
  
  # mageck test zero value
  "mageck_test_remove_zero_threshold": "0" , 
  
  # The name of the file with teh cell line to be used for copy number variation to normalize CNV-biased sgRNA scores prior to gene ranking.
  "cnv_file": "/nexus/posix0/MAGE-flaski/service/projects/data/CRISPR_Screening/CS_main_pipe/CCLE_copynumber_byGene_2013-12-03.txt",
  
  # The name of the cell line to be used for copy number variation to normalize CNV-biased sgRNA scores prior to gene ranking.
  "cnv_line": "U2OS_BONE",
  
  # mageck mle output in project_folder
  "output_mle": "mageck_output/fastq/mle",
  
  # Weight matrix provided in the SSC source for calculating sgRNA efficiencies. If not given SSC will be skipped.
  "efficiency_matrix": "/SSC0.1/matrix/human_CRISPRi_20bp.matrix",
  
  # if mle matrices are being provided specify the folder where they are located
  "mle_matrices":"/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/raw_renamed",
  
  # library.xlx file as generated during preprocessing
  "library_xlsx":"/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/library.xlsx",
  
  # GMT file used for GSEA
  "gmt_file":"/nexus/posix0/MAGE-flaski/service/projects/data/CRISPR_Screening/CS_main_pipe/msigdb.v7.2.symbols.gmt",
  
  # mageck pathway output in project_folder
  "output_pathway":"mageck_output/fastq/pathway",
  
  # mageck plot output in project_folder
  "output_plot":"mageck_output/fastq/plot",
  
  # mageck vispr output in project_folder
  "output_vispr":"mageck_output/fastq/vispr",
  
  # vispr fastq input in project_folder
  "vispr_fastqc":"/nexus/posix0/MAGE-flaski/service/hpc/home/jboucas/nextflow-crispr-data/fastqc_output",
  
  # vispr species
  "vispr_species":"homo_sapiens",
  
  # vispr assembly
  "vispr_assembly":"hg19",
  
  # mageck flute output in project_folder
  "output_flute":"mageck_output/fastq/flute",
  
  # mageck flute organism 
  "mageckflute_organism":"hsa",
  
  # depmap
  "depmap":"True",
  
  # depmap cell line
  "depmap_cell_line":"ACH-000364"
}
```

## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:
```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag> 
```

