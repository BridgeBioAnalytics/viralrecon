process DIVERSITY_STATS {
    label 'process_low'
    cache false

    // conda (params.enable_conda ? 'bioconda::pysam=0.19.1 conda-forge::screed=1.0.5 conda-forge::pandas=1.4.2 conda-forge::click=8.1.2' : null)
    container "821774515222.dkr.ecr.us-west-1.amazonaws.com/nextflow:compute-diversity"

    input:
    path evolved_sequence_counts

    output:
    path 'aa-diversity-stats.csv'       , emit: aa_diversity_stats
    path 'bubble_plot.png'              , emit: bubble_plot
    path 'position_weights_barplot.png' , emit: position_weights_barplot
    path '*position_weights.csv'        , emit: position_weights
    path 'top_kmers_barplot.png'        , emit: top_kmers_barplot
    path 'venn_diagram.png'             , optional:true, emit: venn_diagram
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    """
    echo $PATH
    which -a python
    compute_seq_diversity.py ${evolved_sequence_counts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
