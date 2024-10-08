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

  // when:
  //   ( ! file("${params.project_folder}/samples_renamed.tsv").exists() )

  script:
  """
#!/usr/local/bin/python3
import pandas as pd
import os

print("Starting")
import sys
sys.stdout.flush()

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
  mat_files=os.listdir(matrices)
  mat_files=[s for s in mat_files if ".xlsx" in s]
  print(mat_files)
  for s in mat_files:
      new_name=s.split(".xlsx")[0]
      if not os.path.isfile("${matrices}/mat."+new_name+".tsv"):
          df=pd.read_excel(matrices+"/"+s)
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
if not os.path.isdir( renamed_folder ):
    os.makedirs(renamed_folder)
# if sample names were given without fastq.gz ending. add it
if not 'fastq.gz' in sampleNames.loc[0,1]:
    sampleNames[2] = sampleNames[1] + ['.fastq.gz']
else:
    sampleNames[2] = sampleNames[1]

print(sampleNames)
for index, row in sampleNames.iterrows():
    if not os.path.exists(renamed_folder+"/"+row[2]):
        os.symlink(raw_folder+"/"+row[0], renamed_folder+"/"+row[2])
        print(renamed_folder+"/"+row[2])
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
print("Complete")
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

    mkdir -p ${params.project_folder}/${output_count}

    echo "library: ${library}\noutput_count: ${params.project_folder}/${output_count}"

    mageck count --pdf-report -l ${library} -n ${params.project_folder}/${output_count}/counts --sample-label "\${samples:1}" --fastq \${input_files}

    echo "Mageck count: Done!"
    """
}

process protest {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label
    val paired
    val control
    val treatment
    val control_sgrna
    val control_gene
    val cnv_line
  
  output:
    val "${params.project_folder}/${output_test}", emit: test_output_folder

  script:
    """
    if [ "${params.mageck_test_remove_zero}" == "none" ]
      then
        mageck_test_remove_zero="--remove-zero ${params.mageck_test_remove_zero}"
        mageck_test_remove_zero_threshold=""
    else
        mageck_test_remove_zero="--remove-zero ${params.mageck_test_remove_zero}"
        mageck_test_remove_zero_threshold="--remove-zero-threshold ${params.mageck_test_remove_zero_threshold}"
    fi
    
    if [ "${control}" != "none"  ]
      then 
        control="-c ${control}"
        pdf=""
    else
        pdf="--pdf-report"
    fi
    # not really sure about the pdf variable

    if [ "${cnv_line}" != "" ]
      then
        cnv_norm="--cnv-norm ${params.cnv_file} --cell-line ${cnv_line}"
    elif [ "${params.cnv_line}" != "null" ]  
      then
        cnv_norm="--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
    else 
      cnv_norm=""
    fi

    if [ "${paired}" == paired ]
      then
        paired_testing="--paired"
    else
      paired_testing=""
    fi

    if [ "${control_sgrna}" != "none" ]
      then
        control_sgrna="--control-sgrna ${control_sgrna}"
    else
      control_sgrna=""
    fi

    if [ "${control_gene}" != "none" ]
      then
        control_gene="--control-gene ${control_gene}"
    else
      control_gene=""
    fi

    mageck test \${pdf} --normcounts-to-file \${mageck_test_remove_zero} \${mageck_test_remove_zero_threshold} -k ${params.project_folder}/${params.output_count}/counts.count.txt -t ${treatment} \${control} -n ${params.project_folder}/${params.output_test}/${label} \${cnv_norm} \${paired_testing} \${control_sgrna} \${control_gene}
    
    echo "Mageck test: Done!"
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
    SSC -l ${params.SSC_sgRNA_size} -m ${params.efficiency_matrix} -i ${params.project_folder}/${params.output_mle}/library.txt -o ${params.project_folder}/${params.output_mle}/library.eff.txt  & sleep 60
    echo "SSC: Done!"
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
import pandas as pd
import numpy as np
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
      if cnv_line != "":
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line {cnv_line}"
      elif "${params.cnv_line}" != "null" : 
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
      else:
        cnv_norm=""
    elif "${params.cnv_line}" != "null" : 
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
  if not os.path.exists("${params.project_folder}/${params.output_mle}/mle_matrix") :
    os.makedirs("${params.project_folder}/${params.output_mle}/mle_matrix")

  matfiles=os.listdir(matrices)
  matfiles=[ s for s in matfiles if "mat." in s and ".tsv" in s ]
  for mat in matfiles:
    label=mat.split(".")[1]
    txt=mat.replace(".tsv", ".txt" )
    df=pd.read_csv(f"{matrices}/{txt}",sep=";", header=None)

    control_sgrna=df.loc[ df[0]=="sgrnas", 1].values[0]
    ## if control_sgrna != "none" : 
    if str(control_sgrna) != str(np.nan) :
      control_sgrna=f"--control-sgrna {control_sgrna}"
    else:
      control_sgrna=""

    control_gene=df.loc[ df[0]=="genes", 1].values[0]
    ## if control_gene != "none" :
    if str(control_gene) != str(np.nan) :
      control_gene=f"--control-gene {control_gene}"
    else:
      control_gene=""

    if "cnv_line" in df[0].tolist():
      cnv_line=df.loc[ df[0]=="cnv_line", 1].values[0]
      if cnv_line != "":
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line {cnv_line}"
      elif "${params.cnv_line}" != "null" : 
        cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
      else:
        cnv_norm=""
    elif "${params.cnv_line}" != "null" : 
      cnv_norm=f"--cnv-norm ${params.cnv_file} --cell-line ${params.cnv_line}"
    else:
      cnv_norm=""

    cmd = (
    f"mageck mle "
    f"-k ${params.project_folder}/${params.output_count}/counts.count.txt "
    f"-d '{matrices}/{mat}' "
    f"-n ${params.project_folder}/${params.output_mle}/mle_matrix/{label}_matrix "
    f"{control_sgrna} {control_gene} {cnv_norm} "
    f"--threads=10 ${sgrna_efficiency}"
    )

    print(f"Testing from matrix file {label}")
    print(f"    cmd: {cmd}")
      
    with open(f"${params.project_folder}/${params.output_mle}/mle_matrix/{label}_matrix.mle.sh", "w" ) as f :
      f.write(cmd) 

Path(f"${params.project_folder}/${params.output_mle}/mle.preprocess.done").touch()
print("mageck MLE pre streps - done!")
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
  echo "mageck MLE: Done!"
  """

}

process promle2 {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val sh_script

  output:
    val "${params.project_folder}/${output_test}", emit: test_output_folder

  script:
  """
  bash ${sh_script}
  echo "mageck MLE: Done!"
  """

}

process integration_plot {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  when:
    ( ! file("${params.project_folder}/${params.output_mle}/mle.integration_plot.done").exists()  )

  input:
    val folder
    val folder_status
    val matrices

  script:
    """
#!/usr/local/bin/python3
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from pathlib import Path

project_folder="${params.project_folder}"
mle="${params.output_mle}"
mle_matrix=mle+"/mle_matrix"

### beta from pairwise
def collect_betas(proj):
    proj_betas=pd.DataFrame(columns=["Gene"])
    proj_beta_files=os.listdir(project_folder+"/"+mle)
    proj_beta_files=[s for s in proj_beta_files if "gene_summary.txt" in s and "matrix" not in s]
    proj_beta_files

    for file in proj_beta_files:
        # print(file)
        df=pd.read_csv(project_folder+"/"+mle+"/"+file, sep="\t")
        df_sub=df[["Gene"]+[s for s in df.columns.tolist() if "|beta" in s]]
        cols_to_rename=[s for s in df_sub.columns.tolist() if "|beta" in s][0]
        df_sub=df_sub.rename(columns={cols_to_rename: 'pairwise'+"("+cols_to_rename+")"})
        proj_betas=proj_betas.merge(df_sub, how="right", on="Gene")
        
    return proj_betas

all_betas=collect_betas(project_folder) ### beta from pairwise

raw_pair=pd.read_csv(project_folder + '/samples.tsv', sep=";", header=None)
matrices="${matrices}"

### beta from matrix
def collect_betas_matrix(proj):
    
    proj_betas=pd.DataFrame(columns=["Gene"])

    proj_beta_files=os.listdir(project_folder+"/"+mle_matrix)
    proj_beta_files=[s for s in proj_beta_files if s==f"{new_name}_matrix.gene_summary.txt"][0] 

    df=pd.read_csv(project_folder+"/"+mle_matrix+"/"+proj_beta_files, sep="\t")
    df_sub=df[["Gene"]+[s for s in df.columns.tolist() if "|beta" in s]]
    
    # Rename the columns using the dictionary
    rename_dict = {col: 'matrix' + '('+ col + ')' for col in df_sub.columns if "|beta" in col }
    df_sub = df_sub.rename(columns=rename_dict)
    
    return df_sub

if matrices :
  mat_files=os.listdir(matrices)
  mat_files=[s for s in mat_files if ".tsv" in s]

  for mat_file in mat_files:
    new_name=mat_file.split('mat.')[1].split('.tsv')[0]

    matrix_betas=collect_betas_matrix(project_folder)  ### beta from matrix

    if not os.path.isfile(project_folder+"/"+mle_matrix+"/pairwise_vesus_matrix_correlationPlots"+new_name+".pdf"):

      raw_matrix=pd.read_csv(f"{matrices}/{mat_file}", sep="\t")

      ### go through columns in matrix (from the 3rd column), get sample names that are '1', pick up the relavent rows in pairwise test tsv file
      ### this way, to avoid too many unrelavent pairs of matrix versus pairwise

      matrix_dict={}
      for i in range(2, raw_matrix.shape[1]):
          column_name = raw_matrix.columns[i]
          samples = raw_matrix[raw_matrix.iloc[:, i] == 1].iloc[:, 0]
          matrix_dict[column_name] = samples.tolist()

      pairs_set = set()

      # Iterate over keys in matrix_dict
      for key in matrix_dict:
          # Iterate over the elements in the value list of the current key
          for element in matrix_dict[key]:
              # Check if the element is in any of the specified columns (2 and 3)
              matches = raw_pair[(raw_pair.iloc[:, 2].str.contains(element)) | (raw_pair.iloc[:, 3].str.contains(element))]
              # If matches found, get the corresponding values from the first column
              if not matches.empty:
                  for value in matches.iloc[:, 0]:
                      pairs_set.add((key, value))

      pairs = [list(pair) for pair in pairs_set]

      # label look like ['matrix(Time|beta)', 'pairwise(JEKO_Cyta|beta)'
      def label_pair(pair):
          first_element = f"matrix({pair[0]}|beta)"
          second_element = f"pairwise({pair[1]}|beta)"
          return [first_element, second_element]

      label_pairs = [label_pair(pair) for pair in pairs]

      #plot
      font = {'family' : 'serif',
              'weight' : 'normal',
              'size'   : 16}
          
      pdf = PdfPages(project_folder+"/"+mle_matrix+"/pairwise_vesus_matrix_correlationPlots"+new_name+".pdf")

      w=4
      l=3
      i=1

      fig = plt.figure(figsize=(30,20))

      for pair in label_pairs:
          x=pair[0]
          y=pair[1]
          df_plot=all_betas[['Gene']+[y]].merge(matrix_betas[['Gene']+[x]], how="right", on="Gene")
          
          ax = fig.add_subplot(l,w,i)

          #ax.set_title("%s" %(gene_name), font)
          sns.scatterplot(x=x, y=y, data=df_plot)
          ax.set_xlabel(x, font)
          ax.set_ylabel(y, font)
          ax.tick_params(labelsize=16)
          #matplotlib.rc('axes', linewidth=2)

          i=i+1

          if i == w*l+1:
              plt.tight_layout()
              plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
              plt.show()
              plt.close()
              i=1
              fig = plt.figure(figsize=(30,20))

      if i != 1:
          plt.tight_layout()
          plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
          plt.show()
          plt.close()

      pdf.close()

Path(f"${params.project_folder}/${params.output_mle}/mle.integration_plot.done").touch()
print("mageck MLE integration_plot - done!")
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
print("Merging and Annotation: Done!")
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
# visualize 
dir.create(paste0("${params.project_folder}/${params.output_pathway}/", label), showWarnings = FALSE)
# Negative
sigNeg = subset(fgseaNegScore, padj <= 0.2)
for(i in 1:nrow(sigNeg)){
  P = as.character(sigNeg[i,'pathway'])
  p = plotEnrichment(gmt_file[[P]], genelist_negScore) + labs(title=P)
  ggsave(paste0("neg.", P, '.pdf'), path = paste0("${params.project_folder}/${params.output_pathway}/", label), width = 10, height = 6)
  
  # exctract gene names for flaski
  plot_data = p[['data']]
  a = plot_data[seq(3, nrow(plot_data), 2), 1][[1]]
  # sort genelist
  rnk <- rank(-genelist_negScore)
  ord <- order(rnk)
  statsAdj <- genelist_negScore[ord]
  
  plot_data[seq(3, nrow(plot_data), 2), 'label'] = names(statsAdj)[a]
  write.xlsx(plot_data, paste0("${params.project_folder}/${params.output_pathway}/", label, "/neg.", P, '.xlsx'), row.names = FALSE)
 }
# Positive
sigPos = subset(fgseaPosScore, padj <= 0.2)
for(i in 1:nrow(sigPos)){
  P = as.character(sigPos[i,'pathway'])
  p = plotEnrichment(gmt_file[[P]], genelist_posScore) + labs(title=P)
  ggsave(paste0("pos.", P, '.pdf'), path = paste0("${params.project_folder}/${params.output_pathway}/", label), width = 10, height = 6)

  # exctract gene names for flaski
  plot_data = p[['data']]
  a = plot_data[seq(3, nrow(plot_data), 2), 1][[1]]
  # sort genelist
  rnk <- rank(-genelist_posScore)
  ord <- order(rnk)
  statsAdj <- genelist_posScore[ord]
  plot_data[seq(3, nrow(plot_data), 2), 'label'] = names(statsAdj)[a]
  write.xlsx(plot_data, paste0("${params.project_folder}/${params.output_pathway}/", label, "/pos.", P, '.xlsx'), row.names = FALSE)
}
print("GSEA pathway analysis: Done.")
sessionInfo()
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
  echo "Taring fgsea output: Done."
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
  echo "mageck plot: done."
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
print("vispr: done.")
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
print("vispr-fastqc: done.")
"""
}

process profluterra {
  maxForks 1
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label
    val cell_lines

  when:
    (  ! file("${params.project_folder}/${params.output_test}/MAGeCKFlute_${label}/FluteRRA_${label}.pdf").exists() )

  script:
  """
#!/usr/bin/Rscript
library(MAGeCKFlute)
library(ggplot2)

if ("${cell_lines}" == "none") {
  if ("${params.depmap_cell_line}" != "") {
    cell_lines = "${params.depmap_cell_line}"
  } else {
    cell_lines = NA
  }
} else {
  cell_lines = "${cell_lines}"
}

gdata_file <-"${params.project_folder}/${params.output_test}/${label}.gene_summary.txt"
gdata = ReadRRA(gdata_file)

sdata_file <- "${params.project_folder}/${params.output_test}/${label}.sgrna_summary.txt"
sdata = ReadsgRRA(sdata_file)

if ( "${params.mageckflute_organism}" == "mmu" ) {
  gdata\$HumanGene = TransGeneID(gdata\$id, fromType = "symbol", toType = "symbol", fromOrg = "mmu", toOrg = "hsa")
  idx = duplicated(gdata\$HumanGene)|is.na(gdata\$HumanGene)
  gdata = gdata[!idx, ]
  FluteRRA(gdata, sdata, proj="${label}", organism="${params.mageckflute_organism}", outdir="${params.project_folder}/${params.output_test}/", omitEssential=FALSE, incorporateDepmap = TRUE, cell_lines = cell_lines)
} else {
  FluteRRA("${params.project_folder}/${params.output_test}/${label}.gene_summary.txt", "${params.project_folder}/${params.output_test}/${label}.sgrna_summary.txt", proj="${label}", organism="${params.mageckflute_organism}", outdir="${params.project_folder}/${params.output_test}/", omitEssential=FALSE, incorporateDepmap = TRUE, cell_lines = cell_lines)
}

ResembleDepmap_man <- function(dd, symbol = "id", score = "Score", lineages = "All",
                           method = c("pearson", "spearman", "kendall")[1]){
  dd = dd[!duplicated(dd[, symbol]), ]
  rownames(dd) = dd[, symbol]
  
  ## Load Depmap data
  depmap_rds = "/nexus/posix0/MAGE-flaski/service/databases/MAGeCKFlute/depmap/Depmap_19Q3.rds"
  if(file.exists(depmap_rds)){
    Depmap_19Q3 = readRDS(depmap_rds)
  }else{
    Depmap_19Q3 = t(read.csv("https://ndownloader.figshare.com/files/24613292", header = TRUE,
                             row.names = 1, stringsAsFactors = FALSE, check.names = FALSE))
    rownames(Depmap_19Q3) = gsub(" .*", "", rownames(Depmap_19Q3))
    saveRDS(Depmap_19Q3, depmap_rds)
  }
  meta_rds = "/nexus/posix0/MAGE-flaski/service/databases/MAGeCKFlute/depmap/Depmap_sample_info.rds"
  if(file.exists(meta_rds)){
    sampleinfo = readRDS(meta_rds)
  }else{
    sampleinfo = read.csv("https://ndownloader.figshare.com/files/24613394",
                          row.names = 1, header = TRUE, stringsAsFactors = FALSE)
    saveRDS(sampleinfo, meta_rds)
  }
  if(!"all" %in% tolower(lineages)){
    idx = sampleinfo\$lineage%in%tolower(lineages)
    idx = colnames(Depmap_19Q3)%in%rownames(sampleinfo)[idx]
    if(sum(idx)>5){
      Depmap_19Q3 = Depmap_19Q3[, idx]
    }else{ warning("Less than 5 cell lines are avaible, so ignore lineage setting.")}
  }
  ## Explore the relationship
  genes = intersect(rownames(dd), rownames(Depmap_19Q3))
  if(length(genes)<10) stop("Invalid gene symbols.")
  if(method %in% c("pearson", "spearman", "kendall")){
    similarity = apply(Depmap_19Q3[genes,], 2, function(x){
      tmp = cor.test(x, dd[genes, score], method = method, na.action=na.omit)
      c(tmp\$estimate, tmp\$p.value)
    })
  }else{
    stop("Invalid distance measure!!!")
  }
  similarity = as.data.frame(t(similarity))
  colnames(similarity) = c("estimate", "p.value")
  rownames(similarity) = sampleinfo[colnames(Depmap_19Q3), 1]
  similarity = similarity[order(-similarity\$estimate), ]
  return(similarity)
}

if ( "${params.mageckflute_organism}" == "mmu" ) {
  depmap_sim <- ResembleDepmap_man(gdata, symbol = "HumanGene", score = "Score", lineages = "All", method = c("pearson", "spearman", "kendall")[1])
} else {
  depmap_sim <- ResembleDepmap_man(gdata, symbol = "id", score = "Score", lineages = "All", method = c("pearson", "spearman", "kendall")[1])
}

write.table(depmap_sim, file="${params.project_folder}/${params.output_test}/MAGeCKFlute_${label}/FluteRRA_${label}.depmap.similarity.tsv", sep="\\t", quote = FALSE, row.names = TRUE)

print("FluteRRA: Done.")
  """
}

process proflutemle {
  maxForks 1
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label
    val cell_lines
  when:
    (  ! file("${params.project_folder}/${params.output_mle}/depmap/MAGeCKFlute_${label}/FluteMLE_${label}_cell_cycle.pdf").exists() )

  script:
  """
#!/usr/bin/Rscript
library(MAGeCKFlute)
library(ggplot2)
print("${cell_lines}")

if ("${cell_lines}" == "none") {
  if ("${params.depmap_cell_line}" != "") {
    cell_lines = "${params.depmap_cell_line}"
  } else {
    cell_lines = NA
  }
} else {
  cell_lines = "${cell_lines}"
}

gdata_file <-"${params.project_folder}/${params.output_mle}/${label}.gene_summary.txt"
gdata = ReadBeta(gdata_file)

if ( "${params.mageckflute_organism}" == "mmu" ) {
  gdata\$HumanGene = TransGeneID(gdata\$Gene, fromType = "symbol", toType = "symbol", fromOrg = "mmu", toOrg = "hsa")
  idx = duplicated(gdata\$HumanGene)|is.na(gdata\$HumanGene)
  gdata = gdata[!idx, ]
  FluteMLE(gdata, treatname="${label}", ctrlname="Depmap", proj="${label}", organism="${params.mageckflute_organism}", outdir="${params.project_folder}/${params.output_mle}/depmap", omitEssential = FALSE,  incorporateDepmap=TRUE , cell_lines = cell_lines)
} else {
  FluteMLE("${params.project_folder}/${params.output_mle}/${label}.gene_summary.txt", treatname="${label}", ctrlname="Depmap", proj="${label}", organism="${params.mageckflute_organism}", outdir="${params.project_folder}/${params.output_mle}/depmap", omitEssential = FALSE,  incorporateDepmap=TRUE , cell_lines = cell_lines  )
}

print("FluteMLE: Done.")
"""
}

process promagecku {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label
    val control
    val treatment

  script:
"""
#!/usr/local/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import random
from random import shuffle
import subprocess
import shlex
import sys

def execute(command):
    print(command)
    subprocess.call(shlex.split(command))

def product_threshold_fdr(df,fdr = 0.05):
    maxi = abs(df['product']).max()
    for pro in np.arange(0,maxi,0.1):
        df_thres = df[abs(df['product'])>pro]
        if (1.0 * len(df_thres[df_thres['index'].str.contains('NTC')]) / len(df_thres)) < fdr:
            break
    return pro, df_thres

def rank_test_(gene,ntc_sgRNA_p,df_targeting):
    df_gene = df_targeting[df_targeting['Gene'] == gene].sort_values('p.twosided')
    lfc = df_gene.iloc[:3]['LFC'].mean()
    x, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:3],ntc_sgRNA_p,alternative='two-sided')
    return f"{lfc},{pvalue}" 

def rank_test(df, control_genes ):
    df = df[df['treat_mean'] > 20]

    #df_ntc = df[df['Gene'].str.contains('${params.nontargeting_tag}')]
    #df_targeting =  df[~df['Gene'].str.contains('${params.nontargeting_tag}')]  

    df_ntc = df[df['Gene'].isin(control_genes)]
    df_targeting =  df[~df['Gene'].isin(control_genes)]

    ntc_sgRNA_p = list(df_ntc['p.twosided'])
    ntc_sgRNA_p_lfc = list(zip(list(df_ntc['p.twosided']),list(df_ntc['LFC'])))
    genes = df_targeting['Gene'].unique()
    num_of_genes = len(genes)

    # parallel approach
    gene_lfc_p = {}
    genes=pd.DataFrame({"gene":genes})

    def rank_test_(gene,ntc_sgRNA_p=ntc_sgRNA_p,df_targeting=df_targeting):
        df_gene = df_targeting[df_targeting['Gene'] == gene].sort_values('p.twosided')
        lfc = df_gene.iloc[:3]['LFC'].mean()
        try:
            x, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:3],ntc_sgRNA_p,alternative='two-sided')
        except Exception as error :
            print("##################")
            print("gene",gene, "p values", list(df_gene['p.twosided'])[:3], "ntc_sgrna_p", ntc_sgRNA_p  )
            print("##################")
            x, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:3],ntc_sgRNA_p,alternative='two-sided')
  
        return f"{lfc},{pvalue}" 

    genes["lfc_pvalue"]=genes["gene"].apply(lambda gene : rank_test_(gene) )
    lfc_pvalue=[ [ float( s.split(",")[0] ),   float( s.split(",")[1] ) ] for s in  genes["lfc_pvalue"].tolist() ]
    genes=genes["gene"].tolist()
    # gene_lfc_p=dict(zip(genes, lfc_pvalue ))

    # serial approach
    #for gene in genes:
    #    df_gene = df_targeting[df_targeting['Gene'] == gene].sort_values('p.twosided')
    #    lfc = df_gene.iloc[:3]['LFC'].mean()
    #    print( list(df_gene['p.twosided'])[:3], ntc_sgRNA_p )
    #    x, pvalue = mannwhitneyu(list(df_gene['p.twosided'])[:3],ntc_sgRNA_p,alternative='two-sided')
    #    gene_lfc_p[gene] = [lfc,pvalue]

    # parallel approach
    random.seed(10)
    def ntc_i( ntc_sgRNA_p=ntc_sgRNA_p, ntc_sgRNA_p_lfc=ntc_sgRNA_p_lfc ):
        shuffle(ntc_sgRNA_p_lfc)
        ntc_selected = ntc_sgRNA_p_lfc[:5]
        ntc_selected_p = [i[0] for i in ntc_selected]
        ntc_lfc = np.mean([i[1] for i in sorted(ntc_selected, key=lambda x: x[0])][:3])
        x, ntc_pvalue = mannwhitneyu(ntc_selected_p,ntc_sgRNA_p,alternative='two-sided')
        return f"{ntc_lfc},{ntc_pvalue}" 

    num_of_genes=pd.DataFrame( { "j":range(num_of_genes) } )
    num_of_genes["lfc_pvalue"]=num_of_genes["j"].apply(lambda x: ntc_i() )
    lfc_pvalue_=[ [ float( s.split(",")[0] ),   float( s.split(",")[1] ) ] for s in  num_of_genes["lfc_pvalue"].tolist() ]
    num_of_genes=num_of_genes["j"].tolist()
    num_of_genes=[ 'NTC_' + str(j) for j in num_of_genes ]
    
    gene_lfc_p=dict(zip( genes + num_of_genes , lfc_pvalue + lfc_pvalue_ ))

    # serial approach
    #random.seed(10)
    #for j in range(num_of_genes):
    #    shuffle(ntc_sgRNA_p_lfc)
    #    ntc_selected = ntc_sgRNA_p_lfc[:5]
    #    ntc_selected_p = [i[0] for i in ntc_selected]
    #    ntc_lfc = np.mean([i[1] for i in sorted(ntc_selected, key=lambda x: x[0])][:3])
    #    x, ntc_pvalue = mannwhitneyu(ntc_selected_p,ntc_sgRNA_p,alternative='two-sided')
    #    gene_lfc_p['NTC_' + str(j)] = [ntc_lfc, ntc_pvalue]


    return gene_lfc_p

#print ('fdr: ')
#fdr = float(raw_input('-->'))
fdr = float("${params.magecku_fdr}")

#print ('counts file: ')
#counts_file = raw_input('-->')
counts_file = "${params.project_folder}/${params.output_count}/counts.count.txt"

#print ('control_group label (if more than 1, seperate by ,): ')
#control_group = raw_input('-->')
control_group = "${control}"

#print ('counts threshold_control_groups: ')
#counts_thres_control = float(raw_input('-->'))
counts_thres_control = float("${params.magecku_threshold_control_groups}")

#print ('treatment_group label (if more than 1, seperate by ,): ')
#treatment_group = raw_input('-->')
treatment_group = "${treatment}"

#print ('counts threshold_treatment_groups: ')
#counts_thres_treatment = float(raw_input('-->'))
counts_thres_treatment =  float("${params.magecku_threshold_treatment_groups}")

#print ('output folder: ')
#output_folder = raw_input('-->')
output_folder = "${params.project_folder}/${params.output_magecku}"

#print ('comparison name: ')
#output_name = raw_input('-->')
output_name = "${label}"

print ('path of gene list to label: (enter 0 for not labeling genes) ')
#gene_list_path = raw_input('-->')
gene_list_path = 0

#print ('label all hit genes? (0 or 1) ')
#label_all_sig_genes = raw_input('-->')
label_all_sig_genes = 1 

#print ('plot( 0 or 1 )?: ')
#make_plot = eval(raw_input('-->'))
make_plot = eval("1")

if str(gene_list_path) == '0':
    genes_to_label = []
else:
    with open(gene_list_path.strip(),'r') as f:
        genes_to_label = [i.strip() for i in f.readlines()]

# thresholding
control_groups = str(control_group).split(',')
control_thres = counts_thres_control
treatment_groups = str(treatment_group).split(',')
treatment_thres = counts_thres_treatment
df = pd.read_table(counts_file)
df_thres = df[(df[control_groups] > control_thres).all(axis=1)&(df[treatment_groups] > treatment_thres).all(axis=1)]
df_thres.to_csv(output_folder+'/%s_thresholded_counts.txt'%output_name,sep='\\t',index=False)

# MAGeCK
print("running MAGeCK.....")
execute("mageck test -k " + output_folder+'/%s_thresholded_counts.txt'%output_name + " -t "+ treatment_group + " -c " + control_group + " -n " + output_folder + "/" + output_name + " --pdf-report")
print("running u test.....")
# u test
df_mageck = pd.read_table(output_folder + "/" + output_name + '.sgrna_summary.txt')
df_library=pd.read_excel( "${params.project_folder}/library.xlsx" )
control_genes=df_library[df_library['Annotation'].str.contains('${params.nontargeting_tag}')]["gene_ID"].tolist()
df = pd.DataFrame(rank_test(df_mageck, control_genes)).T
df.columns = ['epsilon','pvalue']
df.reset_index(inplace=True)
df['gene'] = df['index'].apply(lambda x:x.split('_')[0])
df['epsilon'] = -df['epsilon']
df_ntc = df[df['gene']=='NTC']
df['epsilon'] = df['epsilon'] - df_ntc['epsilon'].median()
df['product'] = df['epsilon'] * (-np.log10(df['pvalue']))
df.sort_values('product',ascending=False)


# FDR
thres, df_hits = product_threshold_fdr(df,fdr)
df.sort_values('product',ascending=False).to_csv(output_folder + '/' + output_name + '_all_genes.csv',index=False)
df_hits.sort_values('product',ascending=False).to_csv(output_folder + '/' + output_name + '_fdr%s_product%s_hits.csv'%(fdr,thres),index=False)
df_ntc = df[df['index'].str.contains('NTC')]

########### volcano plot
if make_plot == 1:

    print("plotting....")
    npg = ["#E64B35B2","#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2"]
    plt.figure(figsize=[10,8])
    df_pos = df_hits[df_hits['epsilon']>0]
    df_neg = df_hits[df_hits['epsilon']<0]
    plt.scatter(df_pos['epsilon'],-np.log10(df_pos['pvalue']),c="#DC0000FF",s=5,label='Positive hits')
    plt.scatter(df_neg['epsilon'],-np.log10(df_neg['pvalue']),c="#3C5488FF",s=5,label='Negative hits')
    plt.scatter(df['epsilon'],-np.log10(df['pvalue']),c='#F39B7FFF',s=5,label='Other genes')
    plt.scatter(df_pos['epsilon'],-np.log10(df_pos['pvalue']),c="#DC0000FF",s=5,label=None)
    plt.scatter(df_neg['epsilon'],-np.log10(df_neg['pvalue']),c="#3C5488FF",s=5,label=None)
    plt.scatter(df_ntc['epsilon'],-np.log10(df_ntc['pvalue']),c='grey',s=5,label="Negative control")
    genes=list(df['gene'])
    phenotype = list(df['epsilon'])
    p_value = list(-np.log10(df['pvalue']))
    i = 0
    for x, y, s in zip(phenotype, p_value,genes):
        if s in genes_to_label:
            plt.annotate(s,(x,y),fontsize=16)
            if i == 0:
                plt.scatter(x,y,c='darkgreen',s=20,label = 'Genes of interest')
            else:
                 plt.scatter(x,y,c='darkgreen',s=20)
            i = 1
    if str(label_all_sig_genes) == '1':
        genes=list(df_hits['gene'])
        phenotype = list(df_hits['epsilon'])
        p_value = list(-np.log10(df_hits['pvalue']))
        for x, y, s in zip(phenotype, p_value,genes):
            plt.annotate(s,(x,y),fontsize=16)
        
    x = np.arange(0.01,10,0.01)
    y = [thres / i for i in x]
    plt.plot(x,y,'--', c= 'k',)
    plt.plot(-x,y,'--',c = 'k',label='FDR = %s'%fdr)
    lim = max(abs(df['epsilon'].min() - 1),(df['epsilon'].max() + 2))
    plt.xlim(-lim,lim)
    plt.ylim(0,-np.log10(df['pvalue'].min()) + 0.5)
    plt.legend(loc = 1,fontsize='large',fancybox=True)
    plt.xlabel('Phenotype',fontsize = 14)
    plt.ylabel('-log10 P',fontsize = 14)
    plt.title(output_name,fontsize=18)
    plt.savefig(output_folder + "/" + output_name + '_volcano_plot.pdf')
    plt.show()
print("magecku: done.")
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
  if ( ! file("${params.project_folder}/${params.output_test}").isDirectory() ) {
    file("${params.project_folder}/${params.output_test}").mkdirs()
  }

  rows=Channel.fromPath("${params.samples_tsv}", checkIfExists:true).splitCsv(sep:';')
  rows=rows.filter{ ! file( "${params.project_folder}/${params.output_test}/${it[0]}.gene_summary.txt" ).exists() }
  label=rows.flatMap { n -> n[0] }
  paired=rows.flatMap { n -> n[1] }
  control=rows.flatMap { n -> n[2] }
  control=control.map{ "$it".replace(".fastq.gz","") }
  treatment=rows.flatMap { n -> n[3] }
  treatment=treatment.map{ "$it".replace(".fastq.gz","") }
  control_sgrna=rows.flatMap { n -> n[4] }
  control_gene=rows.flatMap { n -> n[5] }
  cnv_line=rows.flatMap { n -> n[6] }
  cnv_fake=rows.flatMap { n -> "none" }
  cnv_line=cnv_line.concat(cnv_fake)

  protest( label, paired, control, treatment, control_sgrna, control_gene, cnv_line )
  merge_sumaries( "${params.project_folder}/${params.output_test}/", protest.out.collect() )
}

workflow magecku {
  if ( 'output_magecku' in params.keySet() ) {
    if ( ! file("${params.project_folder}/${params.output_magecku}").isDirectory() ) {
      file("${params.project_folder}/${params.output_magecku}").mkdirs()
    }
    rows=Channel.fromPath("${params.samples_tsv}", checkIfExists:true).splitCsv(sep:';')
    rows=rows.filter{ ! file( "${params.project_folder}/${params.output_magecku}/${it[0]}_volcano_plot.pdf" ).exists() }
    label=rows.flatMap { n -> n[0] }
    control=rows.flatMap { n -> n[2] }
    control=control.map{ "$it".replace(".fastq.gz","") }
    treatment=rows.flatMap { n -> n[3] }
    treatment=treatment.map{ "$it".replace(".fastq.gz","") }
    // non-targeting needs to be written on the gene name of the non-targeting sgRNAs as they are used as controls
    promagecku ( label, control, treatment )

  }

}

// workflow mageck_ssc {
//   if ( 'efficiency_matrix' in params.keySet() ) {
//     pressc(  )
//     prossc( pressc.out.collect() )
//   }
// }

workflow mageck_premle {
  if ( "${params.skip_mle}" != "True" ) {

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
}

workflow mageck_mle {

  // rows=Channel.fromPath("${params.samples_tsv}", checkIfExists:true).splitCsv(sep:';')
  // rows=rows.filter{ ! file( "${params.project_folder}/${params.output_mle}/${it[0]}.gene_summary.txt" ).exists() }
  // label=rows.flatMap { n -> n[0] }
  // paired=rows.flatMap { n -> n[1] }
  // control=rows.flatMap { n -> n[2] }
  // control=control.map{ "$it".replace(".fastq.gz","") }
  // treatment=rows.flatMap { n -> n[3] }
  // treatment=treatment.map{ "$it".replace(".fastq.gz","") }
  // control_sgrna=rows.flatMap { n -> n[4] }
  // control_gene=rows.flatMap { n -> n[5] }
  // cnv_line=rows.flatMap { n -> n[6] }
  // cnv_fake=rows.flatMap { n -> "none" }
  // cnv_line=cnv_line.concat(cnv_fake)

  // matrix_control="1,0"
  // matrix_target="1,1"
  // betas="base,"+label}
  // include_samples=control+","+treatment

  // controls=control.replace(",", " ")
  if ( "${params.skip_mle}" != "True" ) {
    if ( 'output_mle' in params.keySet() ) {
      data = channel.fromPath( "${params.project_folder}/${params.output_mle}/*mle.sh" )
      data = data.filter{ ! file( "$it".replace(".mle.sh", ".sgrna_summary.txt") ).exists() }
      promle( data )
          
      if ( 'mle_matrices' in params.keySet() ) {
        data = channel.fromPath( "${params.project_folder}/${params.output_mle}/mle_matrix/*mle.sh" )
        data = data.filter{ ! file( "$it".replace(".mle.sh", ".sgrna_summary.txt") ).exists() }
        promle2( data )
        matrices="${params.mle_matrices}"
        integration_plot("${params.project_folder}/${params.output_mle}/mle_matrix/", promle2.out.collect(), matrices)
      }

      merge_sumaries( "${params.project_folder}/${params.output_mle}/", promle.out.collect() )
    }
  }
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

  rows=Channel.fromPath("${params.samples_tsv}", checkIfExists:true).splitCsv(sep:';')
  depmap_cell_line=rows.flatMap { n -> n[7] }.filter { it != null && it != "" }
  cnv_fake=rows.flatMap { n -> "none" }
  depmap_cell_line=depmap_cell_line.concat(cnv_fake)

  profluterra(labels_test, depmap_cell_line)

  if ( "${params.skip_mle}" != "True" ) {
    if ( 'depmap' in params.keySet()  ) {
      if ( ! file("${params.project_folder}/${params.output_mle}/depmap").isDirectory() ) {
        file("${params.project_folder}/${params.output_mle}/depmap").mkdirs()
      }

      labels_mle=channel.fromPath( "${params.project_folder}/${params.output_mle}/*.gene_summary.txt" )
      labels_mle=labels_mle.map{ "$it.baseName" }
      labels_mle=labels_mle.map{ "$it".replace(".txt","").replace(".gene_summary","") }
      proflutemle(labels_mle, depmap_cell_line)
      
    }
  }
}


