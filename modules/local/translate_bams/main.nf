process TRANSLATE_BAM {
    tag "$fasta"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::pysam=0.19.1 conda-forge::screed=1.0.5 conda-forge::pandas=1.4.2 conda-forge::click=8.1.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    path bams
    path fasta

    output:
    path 'evolved_sequence_counts.csv'      , emit: evolved_7mer_counts
    path 'insertion_summary_percentages.csv', emit: insertion_summary_percentages
    path 'insertion_summary.csv'            , emit: insertion_summary
    path 'cigar_strings.csv'                , emit: cigar_strings
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    """
    translate_bams.py \\
        --fasta ${fasta} \\
        ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
