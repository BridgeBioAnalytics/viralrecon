process TRANSLATE_BAMS {
    label 'process_low'

    // conda (params.enable_conda ? 'bioconda::pysam=0.19.1 conda-forge::screed=1.0.5 conda-forge::pandas=1.4.2 conda-forge::click=8.1.2' : null)
    container "821774515222.dkr.ecr.us-west-1.amazonaws.com/nextflow:translate-bams"

    input:
    path(bams)
    path(bais)
    path fasta
    val directed_evolution_sequence

    output:
    path 'evolved_sequence_counts.csv'                  , emit: evolved_7mer_counts
    path 'evolved_sequence_counts_top_10_per_sample.csv', emit: evolved_7mer_counts_top_10
    path 'insertion_summary.csv'                        , emit: insertion_summary
    path 'cigar_strings.csv'                            , emit: cigar_strings
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    """
    translate_bams.py \\
        --directed-evolution-sequence ${directed_evolution_sequence} \\
        --fasta ${fasta} \\
        ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
