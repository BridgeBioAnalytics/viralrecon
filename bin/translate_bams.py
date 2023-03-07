#!/usr/bin/env python

import os
from asyncio.log import logger
from collections import Counter

import click
import pandas as pd
import pysam
import screed

protein = {
    "TTT": "F",
    "CTT": "L",
    "ATT": "I",
    "GTT": "V",
    "TTC": "F",
    "CTC": "L",
    "ATC": "I",
    "GTC": "V",
    "TTA": "L",
    "CTA": "L",
    "ATA": "I",
    "GTA": "V",
    "TTG": "L",
    "CTG": "L",
    "ATG": "M",
    "GTG": "V",
    "TCT": "S",
    "CCT": "P",
    "ACT": "T",
    "GCT": "A",
    "TCC": "S",
    "CCC": "P",
    "ACC": "T",
    "GCC": "A",
    "TCA": "S",
    "CCA": "P",
    "ACA": "T",
    "GCA": "A",
    "TCG": "S",
    "CCG": "P",
    "ACG": "T",
    "GCG": "A",
    "TAT": "Y",
    "CAT": "H",
    "AAT": "N",
    "GAT": "D",
    "TAC": "Y",
    "CAC": "H",
    "AAC": "N",
    "GAC": "D",
    "TAA": "*",
    "CAA": "Q",
    "AAA": "K",
    "GAA": "E",
    "TAG": "*",
    "CAG": "Q",
    "AAG": "K",
    "GAG": "E",
    "TGT": "C",
    "CGT": "R",
    "AGT": "S",
    "GGT": "G",
    "TGC": "C",
    "CGC": "R",
    "AGC": "S",
    "GGC": "G",
    "TGA": "*",
    "CGA": "R",
    "AGA": "R",
    "GGA": "G",
    "TGG": "W",
    "CGG": "R",
    "AGG": "R",
    "GGG": "G",
    # Degenerate codons
    "GCN": "A",
    "GTN": "V",
    "GGN": "G",
    "TCN": "S",
    "CTN": "L",
    "CCN": "P",
    "CGN": "R",
    "ACN": "T",
}


def translate(dna):
    protein_sequence = ""

    # Generate protein sequence
    range_max = len(dna)
    for i in range(0, len(dna) - (3 + len(dna) % 3) + 1, 3):
        # if protein[dna_sequence[i:i+3]] == "STOP" :
        #     break
        try:
            protein_sequence += protein[dna[i: i + 3]]
        except KeyError:
            raise ValueError(
                "One or more N present in non-degenerate position -- Cannot translate"
            )
    return protein_sequence


class InsertionCounter:
    def __init__(self, bams, directed_evolution_interval, prefix):

        self.bams = bams
        self.directed_evolution_interval = directed_evolution_interval
        self.prefix = prefix
        self.sequence_counter = {}
        self.non_translated_counter = {}
        self.cigar_strings = {}

    def count_27nt_insertions(self):

        for bam in self.bams:
            sample_name = os.path.basename(bam).split(".")[0]
            self.sequence_counter[sample_name] = Counter()
            self.non_translated_counter[sample_name] = Counter()
            self.cigar_strings[sample_name] = Counter()

            read_id_to_classification = {}

            # create a bam files
            bamfile = pysam.AlignmentFile(bam, "rb")

            # For all reads covering the directed evolution interval
            self.classify_reads_in_interval(bamfile, read_id_to_classification, sample_name)

            for read in bamfile:
                if read.query_name not in read_id_to_classification:
                    self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                         category="doesnt_overlap_nns_region")
            read_id_to_classification_series = pd.Series(read_id_to_classification, name=f'{sample_name}_per_read')
            read_id_to_classification_series.to_csv(f'{sample_name}__per_read_categories.csv')

        cigar_strings_df, non_translated_df, sequence_counter_df = self._make_dataframes(sample_name)

        return sequence_counter_df, non_translated_df, cigar_strings_df

    def _make_dataframes(self, sample_name):
        sequence_counter_df = (
            pd.DataFrame(self.sequence_counter).fillna(0).astype(int).sort_index(ascending=True)
        )
        # Sort according to first item, arbitrarily
        sequence_counter_df = sequence_counter_df.sort_values(
            sequence_counter_df.columns[0], ascending=False
        )
        sequence_counter_df.index.names = ["aa_7mer", "nt_21mer"]
        non_translated_df = pd.DataFrame(self.non_translated_counter).fillna(0).astype(int)
        cigar_strings_df = pd.DataFrame(self.cigar_strings).sort_values(
            sample_name, ascending=False
        )
        return cigar_strings_df, non_translated_df, sequence_counter_df

    def classify_reads_in_interval(self, bamfile, read_id_to_classification, sample_name):
        for i, read in enumerate(bamfile.fetch('amplicon', *self.directed_evolution_interval)):
            self.cigar_strings[sample_name][read.cigarstring] += 1

            if read.mapping_quality < 30:
                self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                     category="low_quality_mapping")
                continue

            try:
                interval_positions = [
                    read.get_reference_positions().index(x)
                    for x in self.directed_evolution_interval
                ]
            except ValueError:
                self._check_deletions_wildtype(read, read_id_to_classification, sample_name)

                # Move on to next read
                continue

            if all(x > 1 for x in interval_positions):
                i = interval_positions[0]
                j = interval_positions[1]
                evolved_seq = read.query_alignment_sequence[i:j]
                if len(evolved_seq) != 21:
                    self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                         category=f"insertion_{len(evolved_seq)}nt")
                    continue
                # print(evolved_seq)
                try:
                    evolved_seq_translated = translate(evolved_seq)
                except ValueError:
                    self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                         category="21nt_insertion_too_many_ns")
                    continue
                assert len(evolved_seq_translated) == 7
                pair = (evolved_seq_translated, evolved_seq)
                self.sequence_counter[sample_name][pair] += 1
                self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                     category="21nt_insertion")
            else:
                self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                     category="not_all")
        else:
            self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                 category="low_quality_mapping")

    def _check_deletions_wildtype(self, read, read_id_to_classification, sample_name):
        deletions = [
            (operation, length)
            for (operation, length) in read.cigartuples
            if operation == 2
        ]
        if len(deletions) == 1:
            if deletions[0][1] == 21:
                self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                     category="wild_type")
            elif deletions[0][1] < 21:
                self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                     category="insertion_lt_21")
            else:
                self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                     category="insertion_gt_21")

        else:
            self.increment_counter_classify_read(read, read_id_to_classification, sample_name,
                                                 category="insertion_gt_21")

    def increment_counter_classify_read(self, read, read_id_to_classification, sample_name, category):
        self.non_translated_counter[sample_name][category] += 1
        read_id_to_classification[read.query_name] = category


def get_directed_evolution_interval(fasta, evolution_sequence):
    with screed.open(fasta) as records:
        for record in records:
            seq = record["sequence"]
            start = seq.find(evolution_sequence)
            end = start + len(evolution_sequence)
            # Only look at the first sequence
            break
        logger.warning("Only looking at first sequence in fasta file")
    return start, end


def get_top_n_per_column(sequence_counts_df, n=10):
    dfs = []
    for col in sequence_counts_df.columns:
        df = sequence_counts_df.nlargest(n, col)
        dfs.append(df)
    top_n_per_col = pd.concat(dfs)
    return top_n_per_col


@click.command()
@click.argument("bams", nargs=-1)
@click.option(
    "--fasta",
    type=click.Path(),
    help=(
        "Auto-detect the directed evolution start and end by finding the "
        "sequence provided by --directed-evolution-sequence"
    ),
)
@click.option("--directed-evolution-sequence", default="NNSNNSNNSNNSNNSNNSNNS")
@click.option(
    "--evolved-sequence-counts",
    default="evolved_sequence_counts.csv",
    type=click.Path(),
)
@click.option(
    "--evolved-sequence-counts-top-10",
    default="evolved_sequence_counts_top_10_per_sample.csv",
    type=click.Path(),
)
@click.option("--insertion-summary", default="insertion_summary.csv", type=click.Path())
@click.option(
    "--cigar-strings-csv",
    default="cigar_strings.csv",
    type=click.Path(),
)
@click.option("--prefix", default="")
def main(
    bams,
    fasta,
    directed_evolution_sequence,
    evolved_sequence_counts,
    evolved_sequence_counts_top_10,
    insertion_summary,
    cigar_strings_csv,
    prefix,
):
    if prefix:
        evolved_sequence_counts = f"{prefix}__{evolved_sequence_counts}"
        insertion_summary = f"{prefix}__{insertion_summary}"
        cigar_strings_csv = f"{prefix}__{cigar_strings_csv}"

    directed_evolution_interval = get_directed_evolution_interval(
        fasta, directed_evolution_sequence
    )

    insertion_counter = InsertionCounter(bams, directed_evolution_interval, prefix)

    sequence_counter_df, non_translated_df, cigar_strings_df = insertion_counter.count_27nt_insertions()
    sequence_counter_df.to_csv(evolved_sequence_counts)

    sequence_counter_df_top_n = get_top_n_per_column(sequence_counter_df)
    sequence_counter_df_top_n.to_csv(evolved_sequence_counts_top_10)

    # Write the transpose down so the sample ids are the first column (row names)
    non_translated_df_percentages = 100 * non_translated_df / non_translated_df.sum()
    non_translated_df_percentages.index = 'percent_' + non_translated_df_percentages.index

    insertion_summary_df = pd.concat([non_translated_df, non_translated_df_percentages])
    # Write the transpose down so the sample ids are the first column (row names)
    insertion_summary_df.T.to_csv(insertion_summary)

    cigar_strings_df.to_csv(cigar_strings_csv)


if __name__ == "__main__":
    main()
