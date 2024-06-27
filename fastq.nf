#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq = "/home/alfia/Downloads/fastqc/fastq/*.fastq"
params.qualitycontrol_Result = "/home/alfia/Downloads/fastqc/fastqcresult/"
params.trimmomatic_path = "/usr/share/java/"
params.trimming_Result = "/home/alfia/Downloads/fastqc/trimming/"
params.multiqc_Result = "/home/alfia/Downloads/fastqc/multiqc/"

process QualityControl{

    publishDir("${params.qualitycontrol_Result}", mode: 'copy')
    input:
      path fastq
    
    output:
      path "*"

    script:
    """
    fastqc $fastq
    """

}
process multiqc {
publishDir("${params.multiqc_Result}", mode: 'copy')

input:
    path fastqcresult

output:
    path "*"


script:

"""
multiqc $fastqcresult
"""


}

process trimming {
publishDir("${params.trimming_Result}", mode: 'copy')

input:
 path fastq

output:
 path "*"

script:

"""
java -jar ${params.trimmomatic_path}trimmomatic-0.39.jar SE -phred33 $fastq trimmed_fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
"""

}

workflow {
// Execute quality control
fastq_ch = Channel.fromPath(params.fastq)


QualityControl(fastq_ch)


QualityControl.out.view()

// Execute multiqc
multiqc_ch = Channel.fromPath(params.qualitycontrol_Result)

multiqc(multiqc_ch)
multiqc.out.view()

// Execute trimming
trimming_ch = Channel.fromPath(params.fastq)

trimming(trimming_ch)
trimming.out.view()


}