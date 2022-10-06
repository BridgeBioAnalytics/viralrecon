//
// Translate bam file
//

include { TRANSLATE_BAM } from '../../modules/local/translate_bams/main'
include { DIVERSITY_STATS } from '../../modules/local/diversity_stats/main'


workflow LIBRARY_DIVERSITY {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    fasta         // channel: /path/to/genome.fasta

    main:

    ch_versions = Channel.empty()

    TRANSLATE_BAM(bam.collect(), fasta)

    DIVERSITY_STATS(TRANSLATE_BAM.out.evolved_7mer_counts)

    emit:
    evolved_7mer_counts        = TRANSLATE_BAM.out.evolved_7mer_counts
    evolved_7mer_counts_top_10 = TRANSLATE_BAM.out.evolved_7mer_counts_top_10
    insertion_summary          = TRANSLATE_BAM.out.insertion_summary
    cigar_strings              = TRANSLATE_BAM.out.cigar_strings
    aa_diversity_stats         = DIVERSITY_STATS.out.aa_diversity_stats
    position_weights           = DIVERSITY_STATS.out.position_weights
    versions                   = ch_versions                     // channel: [ versions.yml ]
}
