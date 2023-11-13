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

process count_process {
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

    mageck count --pdf-report  -l ${library} -n ${project_folder}/${output_count}/counts --sample-label "\${samples:1}" --fastq \${input_files}
    """
}

process test_preprocess {
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

if not os.path.exists("${params.project_folder}/${output_test}") :
  os.makedirs("${params.project_folder}/${output_test}")

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
    if "${params.cnv_line}" != "none" : 
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

process test_process {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val sh_script

  script:
    """
    bash ${sh_script}
    """
}

// process upload_paths {
//   stageInMode 'symlink'
//   stageOutMode 'move'

//   input:
//     val fastqc_output

//   script:
//   """
//     cd ${params.project_folder}/${fastqc_output}
//     rm -rf upload.txt
//     for f in \$(ls *.html) ; do echo "fastqc \$(readlink -f \${f})" >>  upload.txt_ ; done
//     uniq upload.txt_ upload.txt 
//     rm upload.txt_
//   """
// }

workflow images {
  get_images()
}

// workflow upload {
//   if ( 'fastqc_output' in params.keySet() ) {
//     fastqc_output="${params.fastqc_output}"
//   } else {
//     fastqc_output="fastqc_output"
//   }
//   upload_paths(fastqc_output)
// }

workflow mageck_count {
    count_process( params.library, params.input_count, params.output_count )
}

workflow mageck_pretest {
    test_preprocess( params.samples_tsv, params.output_test , params.mageck_test_remove_zero, params.mageck_test_remove_zero_threshold )
}

workflow mageck_test {
    data = channel.fromPath( "${params.project_folder}/${params.output_test}/*test.sh" )
    data = data.filter{ ! file("$it".replace(".test.sh", ".gene_summary.txt") ).exists() }
    test_process( data )
}
