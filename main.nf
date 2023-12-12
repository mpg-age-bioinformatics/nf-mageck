#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.containers}" == "singularity" ]] ; 

      then

        cd ${params.image_folder}

        if [[ ! -f mageck-e94a4f7.sif ]] ;
          then
            singularity pull mageck-e94a4f7.sif docker://index.docker.io/mpgagebioinformatics/mageck:e94a4f7
        fi

    fi


    if [[ "${params.containers}" == "docker" ]] ; 

      then

        docker pull mpgagebioinformatics/mageck:e94a4f7

    fi

    """

}

process prorefile{
  stageInMode 'symlink'
  stageOutMode 'move'
  input:
    val matrices

  when:
    ( ! file("${params.project_folder}/samples_renamed.tsv").exists() )

  script:
  """
#!/usr/local/bin/python3
import pandas as pd
import os

matrices="${matrices}"

EXC=pd.ExcelFile("${params.reference_file}", engine="openpyxl")

# process library information
df=EXC.parse("library")
df=df.dropna()
cols=df.columns.tolist()
for c in cols:
    df[c]=df[c].apply( lambda x: str(x).replace(" ","_") )
df.to_excel("${params.project_folder}"+"/library.xlsx", index=None)
# df=df[cols[:3]]
# df.columns=["gene_id","UID", "seq"]
df=df[["gene_ID","UID", "seq"]]
df.columns=["gene_id","UID", "seq"]
df.to_csv("${params.project_folder}"+"/library.tsv", sep="\\t", index=None)
df[["UID", "seq", "gene_id"]].to_csv("${params.project_folder}"+"/library.csv", sep=",", index=None)
dups=df[df.duplicated(subset=["seq"], keep=False)]
# if dups were found report them
if len(dups) > 0 :
    dups.to_excel("${params.project_folder}/library.dups.xlsx", index=None)
df=EXC.parse("sampleNames")
df.to_csv("${params.project_folder}/samplesNames.tsv", sep=";", index=False, header=False)

# process samples
df=EXC.parse("samples")
df.to_csv("${params.project_folder}/samples.tsv", sep=";", index=False, header=False)

# if mle model tables in reference file
sheets=[ s for s in EXC.sheet_names if "mle" in s ]
if matrices :
  for s in sheets:
      new_name=s.split("mle.")[-1].replace(" ","_").replace(".","_")
      if not os.path.isfile("${matrices}/matrix."+new_name+".tsv"):
          df=EXC.parse(s)
          matrix=df[~df["Samples"].isin(["sgrnas","genes"])]
          controls=df[df["Samples"].isin(["sgrnas","genes"])]
          matrix.to_csv("${matrices}/mat."+new_name+".tsv",sep="\\t",index=False)
          controls.to_csv("${matrices}/mat."+new_name+".txt",sep=";",index=False,header=False)

sampleNames = pd.read_csv("${params.project_folder}/samplesNames.tsv", header = None, sep = ";")
samples_file = open('${params.project_folder}/samples.tsv', 'r')
samples = samples_file.read()
samples_file.close()
raw_folder = "${params.raw_fastq}"
renamed_folder = "${params.renamed_fastq}"
# if sample names were given without fastq.gz ending. add it
if not 'fastq.gz' in sampleNames.loc[0,1]:
    sampleNames[2] = sampleNames[1] + ['.fastq.gz']
else:
    sampleNames[2] = sampleNames[1]

for index, row in sampleNames.iterrows():
    if not os.path.exists(renamed_folder+"/"+row[2]):
        os.symlink(raw_folder+"/"+row[0], renamed_folder+"/"+row[2])
    # replace names also in other files, this might be unnesseccary later on
    samples = samples.replace(row[1], row[2])

with open('${params.project_folder}/samples_renamed.tsv', 'w') as O:
    O.write(samples)
#make sure renaming does not happen in the first column 
df_sample = pd.read_csv('${params.project_folder}/samples.tsv', sep=';', header=None)
df_rename = pd.read_csv('${params.project_folder}/samples_renamed.tsv', sep=';', header=None)
# Replace the first column in rename with the orignal
df_rename[0] = df_sample[0]
df_rename.to_csv('${params.project_folder}/samples_renamed.tsv', sep=';', index=False, header=False)
  """

}


process procount {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    val library
    path input_count
    val output_count

  when:
    ( ! file("${params.project_folder}/${output_count}/counts.countsummary.txt").exists() )
  
  script:    
    """
    input_files=""
    samples=""
    cd ${input_count}
    if [[ "${params.mageck_counts_type}" == "bowtie" ]] ; then 
        for f in \$(ls *bam) ; do 
            input_files="\${input_files} \${f}"
            samples="\${samples},\${f%.bam}"
        done
    else
        for f in \$(ls *fastq.gz) ; do 
            input_files="\${input_files} \${f}"
            samples="\${samples},\${f%.fastq.gz}"
        done
    fi

    mkdir -p ${project_folder}/${output_count}

    echo "library: ${library}\noutput_count: ${project_folder}/${output_count}"

    mageck count --pdf-report -l ${library} -n ${project_folder}/${output_count}/counts --sample-label "\${samples:1}" --fastq \${input_files}
    """
}

process pretest {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val samples_tsv
    val output_test
    val mageck_test_remove_zero
    val mageck_test_remove_zero_threshold
  
  when:
    ( ! file("${params.project_folder}/${output_test}/test.preprocess.done").exists() )

  script:
  """
#!/usr/local/bin/python3
import os
from pathlib import Path

if not os.path.exists("${params.project_folder}/${params.output_test}") :
  os.makedirs("${params.project_folder}/${params.output_test}")

mageck_test_remove_zero="--remove-zero ${mageck_test_remove_zero}"
mageck_test_remove_zero_threshold="--remove-zero-threshold ${mageck_test_remove_zero_threshold}"

IDS=""
with open("${samples_tsv}","r") as samples :
  for line in samples:
    line=line.split("\\n")[0]
    l=line.split(";")
    label=l[0]
    paired=l[1]
    control=l[2].replace(".fastq.gz","")
    treatment=l[3].replace(".fastq.gz","")
    control_sgrna=l[4]
    control_gene=l[5]
    if control : control=f"-c {control}"
    if len(l) >= 7:
      cnv_line=l[6]
      if cnv_line != "none":
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line {cnv_line}"
      elif "${params.cnv_line}" != "none" : 
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
      else:
        cnv_norm=""
    elif "${params.cnv_line}" != "none" : 
      cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
    else:
      cnv_norm=""
    if paired : 
      paired_testing="--paired"
    else :
      paired_testing=""
    if control_sgrna != "none" : 
      control_sgrna=f"--control-sgrna {control_sgrna}"
    else:
      control_sgrna=""
    if control_gene != "none" :
      control_gene=f"--control-gene {control_gene}"
    else:
      control_gene=""
    if control : 
      pdf=""
    else :
      pdf="--pdf-report"

    cmd=f"mageck test {pdf} --normcounts-to-file {mageck_test_remove_zero} {mageck_test_remove_zero_threshold} -k ${params.project_folder}/${params.output_count}/counts.count.txt -t {treatment} {control} -n ${params.project_folder}/${output_test}/{label} {cnv_norm} {paired_testing} {control_sgrna} {control_gene}"

    print(f"Testing {label}")
    print(f"    control: {control}")
    print(f"    treatment: {treatment}")

    with open(f"${params.project_folder}/${output_test}/{label}.test.sh", "w" ) as f :
      f.write(cmd)

Path(f"${params.project_folder}/${output_test}/test.preprocess.done").touch()
  """

}

process protest {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val sh_script
  
  output:
    val "${params.project_folder}/${output_test}", emit: test_output_folder

  script:
    """
    bash ${sh_script}
    """
}

process pressc {
  stageInMode 'symlink'
  stageOutMode 'move'

  output:
    val "finished", emit: library_txt

  when:
    ( ! file("${params.project_folder}/${params.output_mle}/library.eff.txt").exists()  )

  script:
    """
#!/usr/local/bin/python3
import os
if not os.path.isdir("${params.project_folder}/${params.output_mle}"):
  os.makedirs("${params.project_folder}/${params.output_mle}")
header=True
with open(f"${params.project_folder}/${params.output_mle}/library.txt", "w") as fout:
  with open("${params.library}", "r" ) as fin:
    for l in fin:
      if header:
        l=l.replace("gene_id", "Gene.Symbol")
        l=l.replace("UID","Row.Names")
        l=l.replace("seq","Sequence")
        header=False
      l=l.split("\\n")[0].split(",")
      l=[ l[1],l[0],l[2] ]
      l="\\t".join(l)
      l=l+"\\n"
      fout.write(l)
    """
}

process prossc {
  stageInMode 'symlink'
  stageOutMode 'move'
  errorStrategy 'ignore'

  input:
    val pre_ssc

  script:
    """
    echo "SSC -l ${params.SSC_sgRNA_size} -m ${params.efficiency_matrix} -i ${params.project_folder}/${params.output_mle}/library.txt -o ${params.project_folder}/${params.output_mle}/library.eff.txt"
    SSC -l ${params.SSC_sgRNA_size} -m ${params.efficiency_matrix} -i ${params.project_folder}/${params.output_mle}/library.txt -o ${params.project_folder}/${params.output_mle}/library.eff.txt
    """
}

process premle {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/${params.output_mle}/mle.preprocess.done").exists()  )

  input:
    val sgrna_efficiency
    val matrices

  output:
    val "finished", emit: premle_output

  script:
    """
#!/usr/local/bin/python3
import os
from pathlib import Path

matrices="${matrices}"

if not os.path.exists("${params.project_folder}/${params.output_mle}") :
  os.makedirs("${params.project_folder}/${params.output_mle}")

IDS=""
with open("${params.samples_tsv}","r") as samples :
  for line in samples:
    line=line.split("\\n")[0]
    l=line.split(";")

    label=l[0]
    paired=l[1]

    control=l[2].replace(".fastq.gz","")
    treatment=l[3].replace(".fastq.gz","")

    control_sgrna=l[4]
    if control_sgrna != "none" : 
      control_sgrna=f"--control-sgrna {control_sgrna}"
    else:
      control_sgrna=""

    control_gene=l[5]
    if control_gene != "none" :
      control_gene=f"--control-gene {control_gene}"
    else:
      control_gene=""

    if len(l) >= 7:
      cnv_line=l[6]
      if cnv_line != "none":
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line {cnv_line}"
      elif "${params.cnv_line}" != "none" : 
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
      else:
        cnv_norm=""
    elif "${params.cnv_line}" != "none" : 
      cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
    else:
      cnv_norm=""

    controls=control.split(",")
    treatments=treatment.split(",")

    matrix_control="1,0"
    matrix_target="1,1"
    betas=f"base,{label}"
    include_samples=f"{control},{treatment}"
    print(l[2],l[3])
    print(control)
    print(treatment)
    print(include_samples)

    design_matrix_controls=[ matrix_control for s in controls ]
    design_matrix_treatments=[ matrix_target for s in treatments ]
    design_matrix=design_matrix_controls+design_matrix_treatments
    design_matrix=";".join(design_matrix)

    cmd=f"mageck mle -k ${params.project_folder}/${params.output_count}/counts.count.txt -d '{design_matrix}' -n ${params.project_folder}/${params.output_mle}/{label} -b {betas} -i {include_samples} {control_sgrna} {control_gene} {cnv_norm} --threads=10 ${sgrna_efficiency}"

    print(f"Testing {label}")
    print(f"    control: {control}")
    print(f"    treatment: {treatment}")
    print(f"    cmd: {cmd}")

    with open(f"${params.project_folder}/${params.output_mle}/{label}.mle.sh", "w" ) as f :
      f.write(cmd)

if matrices:
  matfiles=os.listdir(matrices)
  matfiles=[ s for s in matfiles if "mat." in s ]
  for mat in matfiles:
    label=mat.split(".")[1]
    txt=mat.replace(".tsv", ".txt" )
    df=pd.read_csv(f"{matrices}/{txt}",sep=";", header=None)

    control_sgrna=df.loc[ df[0]=="sgrnas", 1].values[0]
    if control_sgrna != "none" : 
      control_sgrna=f"--control-sgrna {control_sgrna}"
    else:
      control_sgrna=""

    control_gene=df.loc[ df[0]=="genes", 1].values[0]
    if control_gene != "none" :
      control_gene=f"--control-gene {control_gene}"
    else:
      control_gene=""

    if "cnv_line" in df[0].tolist():
      cnv_line=df.loc[ df[0]=="cnv_line", 1].values[0]
      if cnv_line != "none":
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line {cnv_line}"
      elif "${params.cnv_line}" != "none" : 
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
      else:
        cnv_norm=""
    elif "${params.cnv_line}" != "none" : 
      cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
    else:
      cnv_norm=""

    cmd=f"mageck mle -k ${params.project_folder}/${params.output_count}/counts.count.txt -d f"'{matrices}/{mat}'" -n ${params.project_folder}/${params.output_mle}/{label}_matrix {control_sgrna} {control_gene} {cnv_norm} --threads=10 ${sgrna_efficiency}"

    print(f"Testing from matrix file {label}")
    print(f"    cmd: {cmd}")
      
    with open(f"${params.project_folder}/${params.output_mle}/{label}_matrix.mle.sh", "w" ) as f :
      f.write(cmd) 

Path(f"${params.project_folder}/${params.output_mle}/mle.preprocess.done").touch()

    """ 
}

process promle {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val sh_script

  output:
    val "${params.project_folder}/${output_test}", emit: test_output_folder

  script:
  """
  bash ${sh_script}
  """

}

process merge_sumaries {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val folder
    val folder_status

  script:
    """
#!/usr/local/bin/python3
import pandas as pd
import os
# get all gene and sgrna files
infolder = "${folder}"
files = os.listdir(infolder)
gene_files = [file for file in files if 'gene_summary.txt' in file]
sgrna_files = [file for file in files if 'sgrna_summary.txt' in file]
# get annotation file
annotation = pd.read_excel("${params.library_xlsx}", engine="openpyxl")
annotation_gene = annotation[['gene_ID', 'Annotation']]
annotation_sgrna = annotation[['UID', 'Annotation']]
annotation_gene.drop_duplicates(inplace = True)
annotation_sgrna.drop_duplicates(inplace = True)
# annotating and merging gene files
genes = dict()
for gene in gene_files:
    print(gene.replace('.gene_summary.txt', ''))
    tmp = pd.read_csv(infolder + gene, sep = "\\t")
    tmp['analysis'] = gene.replace('.gene_summary.txt', '')
    if "id" in tmp.columns:
        right_on = "id"
    else:
        right_on = "Gene"
    tmp = pd.merge(annotation_gene, tmp, left_on = "gene_ID", right_on = right_on, how = 'right')
    tmp.drop(columns=['gene_ID'], inplace = True)
    tmp.to_excel(infolder + gene.replace('.txt', '.annotated.xlsx'), index = False)
    genes[gene.replace('.gene_summary.txt', '')] = tmp
genes_merged = pd.concat(genes)
genes_merged.reset_index(drop=True, inplace = True)
genes_merged.sort_values(by = [right_on, 'analysis'], inplace = True)
# annotating and merging sgrnas
sgrnas = dict()
for sgrna in sgrna_files:
    print(sgrna.replace('.sgrna_summary.txt', ''))
    tmp = pd.read_csv(infolder + sgrna, sep = "\\t")
    tmp['analysis'] = sgrna.replace('.sgrna_summary.txt', '')
    if "sgrna" in tmp.columns:
        right_on = 'sgrna'
    else:
        right_on = "sgRNA"
    tmp = pd.merge(annotation_sgrna, tmp, left_on = "UID", right_on = right_on, how = 'right')
    tmp.drop(columns=['UID'], inplace = True)
    tmp.to_excel(infolder + sgrna.replace('.txt', '.annotated.xlsx'), index = False)
    sgrnas[sgrna.replace('.sgrna_summary.txt', '')] = tmp
sgrnas_merged = pd.concat(sgrnas)
sgrnas_merged.reset_index(drop=True, inplace = True)
sgrnas_merged.sort_values(by = [right_on, 'analysis'], inplace = True)
# Saving merged files
sgrnas_merged.to_csv(infolder +'merged.sgrnas.summary.tsv', sep='\\t', index=False)
genes_merged.to_csv(infolder + 'merged.genes.summary.tsv', sep='\\t', index = False)
sgrnas_merged.to_excel(infolder + 'merged.sgrnas.summary.xlsx', index = False)
genes_merged.to_excel(infolder + 'merged.genes.summary.xlsx', index = False)
    """

}

process propathway {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val gene_ranking_file 

  output:
    val "${gene_ranking_file}"

  script:
  """
#!/usr/bin/Rscript
library(fgsea)
library(data.table)
library(ggplot2)
library(openxlsx)
# input data
label=strsplit("${gene_ranking_file}", "${params.project_folder}/${params.output_test}/")
label=unlist(label)[2]
label=strsplit(label, ".gene_summary.txt")
label=unlist(label)[1]
gene_rank_file = read.delim("${gene_ranking_file}", as.is = TRUE)
gmt_file <- gmtPathways("${params.gmt_file}")
# using score and converting using z-score normalization
# based on negative score
genelist_negScore = gene_rank_file[, 'neg.score']
names(genelist_negScore) = gene_rank_file[,'id']
# z-score
genelist_negScore = (genelist_negScore - mean(genelist_negScore))/sd(genelist_negScore)
#based on positive score
genelist_posScore = gene_rank_file[, 'pos.score']
names(genelist_posScore) = gene_rank_file[, 'id']
# zscore
genelist_posScore = (genelist_posScore - mean(genelist_posScore))/sd(genelist_posScore)
genelist_posScore = sort(genelist_posScore)
# run fgsea
fgseaNegScore <- fgsea(pathways = gmt_file, stats = genelist_negScore, scoreType = 'std', minSize = 5, maxSize = 2000)
fgseaPosScore <- fgsea(pathways = gmt_file, stats = genelist_posScore, scoreType = 'std', minSize = 5, maxSize = 2000)
# add number of genes enriched in pathway
fgseaNegScore[, 'neg.topGenes'] = unlist(lapply(fgseaNegScore[['leadingEdge']], function(x) length(unlist(x))))
fgseaPosScore[, 'pos.topGenes'] = unlist(lapply(fgseaPosScore[['leadingEdge']], function(x) length(unlist(x))))
# add total number of genes in pathway
x = lapply(gmt_file, length)
pathway_genes = as.data.frame(unlist(x))
names(pathway_genes) = 'full_pathway'
pathway_genes[, 'pathway'] = row.names(pathway_genes)
fgseaNegScore = merge(fgseaNegScore, pathway_genes, by = "pathway", all.x = TRUE )
fgseaPosScore = merge(fgseaPosScore, pathway_genes, by = "pathway", all.x = TRUE )
# sort pathways according to adjusted pvalue
fgseaNegScore = fgseaNegScore[order(fgseaNegScore[, 'padj'], decreasing = FALSE),]
fgseaPosScore = fgseaPosScore[order(fgseaPosScore[, 'padj'], decreasing = FALSE),]
# save data
write.xlsx(fgseaNegScore, file=paste0("${params.project_folder}/${params.output_pathway}/", label, ".neg.fgsea.xlsx"), row.names = FALSE)
write.xlsx(fgseaPosScore, file=paste0("${params.project_folder}/${params.output_pathway}/", label, ".pos.fgsea.xlsx"), row.names = FALSE)
# Vizualise 
dir.create(paste0("${params.project_folder}/${params.output_pathway}/", label), showWarnings = FALSE)
# Negative
sigNeg = subset(fgseaNegScore, padj <= 0.2)
#for(i in 1:nrow(sigNeg)){
for(i in 1:2){
  P = as.character(sigNeg[i,'pathway'])
  p = plotEnrichment(gmt_file[[P]], genelist_negScore) + labs(title=P)
  ggsave(paste0("neg.", P, '.pdf'), path = paste0("${params.project_folder}/${params.output_pathway}/", label), width = 10, height = 6)
  
  # exctract gene names for flaski
  #plot_data = p[['data']]
  #a = plot_data[seq(3, nrow(plot_data), 2), 'x']
  ## sort genelist
  #rnk <- rank(-genelist_negScore)
  #ord <- order(rnk)
  #statsAdj <- genelist_negScore[ord]
  
  #plot_data[seq(3, nrow(plot_data), 2), 'label'] = names(statsAdj)[a]
  #write.xlsx(plot_data, paste0("${params.project_folder}/${params.output_pathway}/", label, "/neg.", P, '.xlsx'), row.names = FALSE)
 }
# Positive
sigPos = subset(fgseaPosScore, padj <= 0.2)
#for(i in 1:nrow(sigPos)){
for(i in 1:2){
  P = as.character(sigPos[i,'pathway'])
  p = plotEnrichment(gmt_file[[P]], genelist_posScore) + labs(title=P)
  ggsave(paste0("pos.", P, '.pdf'), path = paste0("${params.project_folder}/${params.output_pathway}/", label), width = 10, height = 6)

  ## exctract gene names for flaski
  #plot_data = p[['data']]
  #a = plot_data[seq(3, nrow(plot_data), 2), 'x']
  ## sort genelist
  #rnk <- rank(-genelist_posScore)
  #ord <- order(rnk)
  #statsAdj <- genelist_posScore[ord]
  #plot_data[seq(3, nrow(plot_data), 2), 'label'] = names(statsAdj)[a]
  #write.xlsx(plot_data, paste0("${params.project_folder}/${params.output_pathway}/", label, "/pos.", P, '.xlsx'), row.names = FALSE)
}

  """

}

process propathtargz {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val folder
    val pathway_status
  
  script:
  """
  cd ${params.project_folder}/${params.output_pathway}/ && \
  find ${folder} -type f -name '*' > temp.${folder}.txt && \
  tar -T temp.${folder}.txt -cvzf ${folder}.tar.gz --remove-files && \
  rm temp.${folder}.txt && \
  rm -rf ${folder}
  """
}

process proplot {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val gene_ranking_file
    val label
  
  script:
  """
  mageck plot -k ${params.project_folder}/${params.output_count}/counts.count.txt -g ${gene_ranking_file} -n ${params.project_folder}/${params.output_plot}/${label}
  """
}

process provispr {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val gene_ranking_file
    val test_type

  script:
"""
#!/usr/local/bin/python3
import os
import shutil
import zipfile

f=os.path.basename("${gene_ranking_file}")
e=f.split('.gene_summary.txt')[0]
if not os.path.isdir("${params.project_folder}/${params.output_vispr}/${test_type}/") :
  os.makedirs("${params.project_folder}/${params.output_vispr}/${test_type}/")
shutil.copy("${gene_ranking_file}", f"${params.project_folder}/${params.output_vispr}/${test_type}/{f}") 

yaml=f'''experiment: {e}.${test_type}\\n\
species: ${params.vispr_species}\\n\
assembly: ${params.vispr_assembly}\\n\
targets:\\n\
    results: /vispr/${test_type}/{f}\\n\
    genes: true\\n\
sgrnas:\\n\
    counts: /vispr/counts.count_normalized.txt\\n\
    mapstats: /vispr/counts.countsummary.txt\\n\
fastqc:\\n\
'''

fastqc_files=os.listdir("${params.vispr_fastqc}") 
fastqc_files=[ s for s in fastqc_files if ".zip" in s ]
for zip_file in fastqc_files:
  sample=zip_file.split("_fastqc.zip")[0]
  folder=zip_file.split(".zip")[0]
  yaml=f'''{yaml}    {sample}:\\n    - /vispr/fastqc/{folder}/fastqc_data.txt\\n'''

yaml_file=f"${params.project_folder}/${params.output_vispr}/{e}.${test_type}.yaml"
with open(yaml_file, "w") as fout:
  fout.write(yaml)
"""
}

process provispr_fastqc {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
"""
#!/usr/local/bin/python3
import os
import shutil
import zipfile

fastqc_files=os.listdir("${params.vispr_fastqc}") 
fastqc_files=[ s for s in fastqc_files if ".zip" in s ]
for zip_file in fastqc_files:
  sample=zip_file.split("_fastqc.zip")[0]
  folder=zip_file.split(".zip")[0]
  if not os.path.isdir("${params.project_folder}/${params.output_vispr}/fastqc/{folder}/") :
    with zipfile.ZipFile(f"${params.vispr_fastqc}/{zip_file}") as myzip:
      myzip.extract(f'{folder}/fastqc_data.txt',path=f"${params.project_folder}/${params.output_vispr}/fastqc/")
"""
}

process profluterra {
   stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label

  script:
  """
#!/usr/local/bin/Rscript
library(MAGeCKFlute)
library(ggplot2)
FluteRRA("${params.project_folder}/${params.output_test}/${label}.gene_summary.txt", "${params.project_folder}/${params.output_test}/${label}.sgrna_summary.txt", proj="${label}", organism="${params.mageckflute_organism}", outdir="${params.project_folder}/${params.output_test}/" )
  """
}

process proflutemle {
   stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label
    val cell_lines

  script:
  """
#!/usr/local/bin/Rscript
library(MAGeCKFlute)
library(ggplot2)
FluteMLE("${params.project_folder}/${params.output_mle}/${label}.gene_summary.txt", treatname="${label}", ctrlname="Depmap", proj="${label}", organism="${params.mageckflute_organism}", outdir="${params.project_folder}/${params.output_mle}/depmap", incorporateDepmap=TRUE ${cell_lines}  )
  """
}


workflow images {
  get_images()
}

workflow pre_process{
   if ( 'mle_matrices' in params.keySet() ) {
    matrices="${params.mle_matrices}"
  } else {
    matrices=""
  }
  if ( ! file("${params.renamed_fastq}").isDirectory() ) {
    file("${params.renamed_fastq}").mkdirs()
  }
  prorefile(matrices)
}

workflow mageck_count {
  procount( params.library, params.input_count, params.output_count )
}

workflow mageck_pretest {
  pretest( params.samples_tsv, params.output_test , params.mageck_test_remove_zero, params.mageck_test_remove_zero_threshold )
}

workflow mageck_test {
  data = channel.fromPath( "${params.project_folder}/${params.output_test}/*test.sh" )
  data = data.filter{ ! file("$it".replace(".test.sh", ".gene_summary.txt") ).exists() }
  protest( data )
  merge_sumaries( "${params.project_folder}/${params.output_test}/", protest.out.collect() )
}

// workflow mageck_ssc {
//   if ( 'efficiency_matrix' in params.keySet() ) {
//     pressc(  )
//     prossc( pressc.out.collect() )
//   }
// }

workflow mageck_premle {
  if ( 'efficiency_matrix' in params.keySet() ) {
    pressc(  )
    prossc( pressc.out.collect() )
    sgrna_efficiency="--sgrna-efficiency ${params.project_folder}/${params.output_mle}/library.eff.txt --sgrna-eff-name-column 1 --sgrna-eff-score-column 3"
  } else {
    sgrna_efficiency=""
  }

  if ( 'mle_matrices' in params.keySet() ) {
    matrices="${params.mle_matrices}"
  } else {
    matrices=""
  }

  premle(sgrna_efficiency, matrices)

}

workflow mageck_mle {
  data = channel.fromPath( "${params.project_folder}/${params.output_mle}/*mle.sh" )
  data = data.filter{ ! file( "$it".replace(".mle.sh", ".sgrna_summary.txt") ).exists() }
  promle( data )
  merge_sumaries( "${params.project_folder}/${params.output_mle}/", promle.out.collect() )
}

workflow mageck_pathway {
  if ( 'gmt_file' in params.keySet() ) {
    if ( ! file("${params.project_folder}/${params.output_pathway}").isDirectory() ) {
      file("${params.project_folder}/${params.output_pathway}").mkdirs()
    }
    data = channel.fromPath( "${params.project_folder}/${params.output_test}/*.gene_summary.txt" )
    data = data.filter{ ! file( "$it".replace("${params.output_test}", "${params.output_pathway}").replace(".gene_summary.txt",".tar.gz") ).exists() }
    // println "${params.project_folder}/${params.output_test}/*.gene_summary.txt"
    // data.subscribe { println "value: $it" }
    propathway( data )
    folders=data.map{ "$it.baseName" }
    folders=folders.map{ "$it".replace(".txt","").replace(".gene_summary","") }
    propathtargz ( folders, propathway.out.collect() )
  }
}

workflow mageck_plot {
  if ( ! file("${params.project_folder}/${params.output_plot}").isDirectory() ) {
    file("${params.project_folder}/${params.output_plot}").mkdirs()
  }
  data = channel.fromPath( "${params.project_folder}/${params.output_test}/*.gene_summary.txt" )
  data = data.filter{ ! file( "$it".replace("${params.output_test}", "${params.output_plot}").replace(".gene_summary.txt",".pdf") ).exists() }
  labels=data.map{ "$it.baseName" }
  labels=labels.map{ "$it".replace(".txt","").replace(".gene_summary","") }
  proplot(data,labels)
}



workflow mageck_vispr {
  if ( ! file("${params.project_folder}/${params.output_vispr}/test").isDirectory() ) {
    file("${params.project_folder}/${params.output_vispr}/test").mkdirs()
  }
  if ( ! file("${params.project_folder}/${params.output_vispr}/mle").isDirectory() ) {
    file("${params.project_folder}/${params.output_vispr}/mle").mkdirs()
  }

  if ( ! file("${params.project_folder}/${params.output_vispr}/counts.count_normalized.txt").exists() ) {
    file("${params.project_folder}/${params.output_count}/counts.count_normalized.txt").copyTo("${params.project_folder}/${params.output_vispr}/counts.count_normalized.txt")
  }
  if ( ! file("${params.project_folder}/${params.output_vispr}/counts.countsummary.txt").exists() ) {
    file("${params.project_folder}/${params.output_count}/counts.countsummary.txt").copyTo("${params.project_folder}/${params.output_vispr}/counts.countsummary.txt")
  }

  summaries = channel.fromPath( "${params.project_folder}/${params.output_test}/*.gene_summary.txt" )
  summaries = summaries.filter{ ! file( "$it".replace("${params.output_test}", "${params.output_vispr}").replace(".gene_summary.txt",".test.yaml") ).exists() }
  summaries_types=summaries.flatMap{ n -> "test"}
  summaries_ = channel.fromPath( "${params.project_folder}/${params.output_mle}/*.gene_summary.txt" )
  summaries_ = summaries_.filter{ ! file( "$it".replace("${params.output_mle}", "${params.output_vispr}").replace(".gene_summary.txt",".mle.yaml") ).exists() }
  summaries_types_=summaries_.flatMap{ n -> "mle"}

  summaries=summaries.concat(summaries_)
  summaries_types=summaries_types.concat(summaries_types_)

  provispr(summaries,summaries_types)

  if ( ! file("${params.project_folder}/${params.output_vispr}/fastqc").isDirectory() ) {
    file("${params.project_folder}/${params.output_vispr}/fastqc").mkdirs()
    provispr_fastqc()
  }

}

workflow mageck_flute {
  labels_test=channel.fromPath( "${params.project_folder}/${params.output_test}/*.gene_summary.txt" )
  labels_test=labels_test.map{ "$it.baseName" }
  labels_test=labels_test.map{ "$it".replace(".txt","").replace(".gene_summary","") }
  labels_test.view()
  profluterra("AJ02_Lib1_Lib2")

  if ( 'depmap' in params.keySet()  ) {
    if ( ! file("${params.project_folder}/${params.output_mle}/depmap").isDirectory() ) {
      file("${params.project_folder}/${params.output_mle}/depmap").mkdirs()

    if ( 'depmap_cell_line' in params.keySet()  ) {    
      depmap_cell_line="${params.depmap_cell_line}".replace(' ', '","')
      depmap_cell_line=', cell_lines=c("'+depmap_cell_line +'")'
    } else {
      depmap_cell_line=', cell_lines = rownames(depmap_similarity)[1], lineages = "All"'
    }

    labels_mle=channel.fromPath( "${params.project_folder}/${params.output_mle}/*.gene_summary.txt" )
    labels_mle=labels_mle.map{ "$it.baseName" }
    labels_mle=labels_mle.map{ "$it".replace(".txt","").replace(".gene_summary","") }

    proflutemle(labels_mle, depmap_cell_line)

    }


  }

}


