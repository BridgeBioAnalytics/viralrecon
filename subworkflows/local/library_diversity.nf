//
// Translate bam file
//

include { TRANSLATE_BAM } from '../../modules/local/translate_bams/main'
include { COMPUTE_DIVERSITY } from '../../modules/local/compute_diversity/main'


workflow LIBRARY_DIVERSITY {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    fasta         // channel: /path/to/genome.fasta
    sizes         // channel: /path/to/genome.sizes
    gff           // channel: /path/to/genome.gff
    bed           // channel: /path/to/primers.bed
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config

    main:

    ch_versions = Channel.empty()

    TRANSLATE_BAM(bam.collect(), fasta)

    COMPUTE_DIVERSITY(TRANSLATE_BAM.out.evolved_7mer_counts)

    emit:

    versions        = ch_versions                     // channel: [ versions.yml ]
}
