{//
// Translate bam file
//

include { TRANSLATE_BAMS } from '../../modules/local/translate_bams/main'
include { DIVERSITY_STATS } from '../../modules/local/diversity_stats/main'


workflow LIBRARY_DIVERSITY {
    take:
    bam                         // channel: [ val(meta), [ bam ] ]
    fasta                       // channel: /path/to/genome.fasta
    directed_evolution_sequence // val("NNSNNSNNSNNSNNSNNSNNS") for example

    main:

    ch_versions = Channel.empty()

    bams = bam.map{it -> it[1]}.collect().dump(tag: "bams_for_translate")
    bais = bam.map{it -> it[2]}.collect().dump(tag: "bais_for_translate")
    TRANSLATE_BAMS(bams, bais, fasta, directed_evolution_sequence)
    ch_versions = ch_versions.mix(TRANSLATE_BAMS.out.versions.first().ifEmpty(null))

    DIVERSITY_STATS(TRANSLATE_BAMS.out.evolved_7mer_counts)
    ch_versions = ch_versions.mix(DIVERSITY_STATS.out.versions.first().ifEmpty(null))

    emit:
    evolved_7mer_counts        = TRANSLATE_BAMS.out.evolved_7mer_counts
    evolved_7mer_counts_top_10 = TRANSLATE_BAMS.out.evolved_7mer_counts_top_10
    insertion_summary          = TRANSLATE_BAMS.out.insertion_summary
    cigar_strings              = TRANSLATE_BAMS.out.cigar_strings
    aa_diversity_stats         = DIVERSITY_STATS.out.aa_diversity_stats
    position_weights           = DIVERSITY_STATS.out.position_weights
    bubble_plot                = DIVERSITY_STATS.out.bubble_plot
    top_kmers_barplot          = DIVERSITY_STATS.out.top_kmers_barplot
    position_weights_barplot   = DIVERSITY_STATS.out.position_weights_barplot
    versions                   = ch_versions                     // channel: [ versions.yml ]
}
}
