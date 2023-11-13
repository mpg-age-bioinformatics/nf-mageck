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
            samples="\${samples},\${f%.fastq}"
        done
    fi

    mkdir -p ${project_folder}/${output_count}

    echo "library: ${library}\noutput_count: ${project_folder}/${output_count}"

    mageck count --pdf-report  -l ${library} -n ${project_folder}/${output_count}/counts --sample-label "\${samples:1}" --fastq \${input_files}
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