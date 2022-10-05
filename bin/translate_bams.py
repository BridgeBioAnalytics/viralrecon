from collections import Counter
import os

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
            protein_sequence += protein[dna[i : i + 3]]
        except KeyError:
            raise ValueError(
                "One or more N present in non-degenerate position -- Cannot translate"
            )
    return protein_sequence


def count_21nt_7aa(bams, directed_evolution_interval):
    sequence_counter = {}
    non_translated_counter = {}
    cigar_strings = {}

    for bam in bams:
        sample_name = os.path.basename(bam).split(".")[0]
        sequence_counter[sample_name] = Counter()
        non_translated_counter[sample_name] = Counter()
        cigar_strings[sample_name] = Counter()

        # create a bam files
        bamfile = pysam.AlignmentFile(bam, "rb")

        for i, read in enumerate(bamfile.fetch()):
            cigar_strings[sample_name][read.cigarstring] += 1

            if read.mapping_quality < 40:
                continue

            try:
                interval_positions = [
                    read.get_reference_positions().index(x)
                    for x in directed_evolution_interval
                ]
            except ValueError:
                deletions = [
                    (operation, length)
                    for (operation, length) in read.cigartuples
                    if operation == 2
                ]
                # print(deletions)
                if len(deletions) == 1:
                    if deletions[0][1] == 21:
                        non_translated_counter[sample_name]["wild_type"] += 1
                    elif deletions[0][1] < 21:
                        non_translated_counter[sample_name]["other"] += 1
                    else:
                        non_translated_counter[sample_name]["other"] += 1
                else:
                    non_translated_counter[sample_name]["other"] += 1

                # Move on to next read
                continue

            if all(x > 1 for x in interval_positions):
                i = interval_positions[0]
                j = interval_positions[1]
                evolved_seq = read.query_alignment_sequence[i:j]
                if len(evolved_seq) != 21:
                    non_translated_counter[sample_name][
                        f"insertion_{len(evolved_seq)}nt"
                    ] += 1
                    continue
                # print(evolved_seq)
                try:
                    evolved_seq_translated = translate(evolved_seq)
                except ValueError:
                    non_translated_counter[sample_name][
                        "21nt_insertion_too_many_ns"
                    ] += 1
                    continue
                assert len(evolved_seq_translated) == 7
                if read.query_name == "A01940:28:GW220807000:1:2101:18177:1063":
                    import pdb

                    pdb.set_trace()
                pair = (evolved_seq_translated, evolved_seq)
                sequence_counter[sample_name][pair] += 1
                non_translated_counter[sample_name]["21nt_insertion"] += 1
        else:
            non_translated_counter[sample_name]["low_quality_mapping"] += 1

    sequence_counter_df = (
        pd.DataFrame(sequence_counter).fillna(0).astype(int).sort_index(ascending=True)
    )
    # Sort according to first item, arbitrarily
    sequence_counter_df = sequence_counter_df.sort_values(
        sequence_counter_df.columns[0], ascending=False
    )
    sequence_counter_df.index.names = ["aa_7mer", "nt_21mer"]

    non_translated_df = pd.DataFrame(non_translated_counter).fillna(0).astype(int)

    cigar_strings_df = pd.DataFrame(cigar_strings).sort_values(
        sample_name, ascending=False
    )

    return sequence_counter_df, non_translated_df, cigar_strings_df


def get_directed_evolution_interval(fasta, evolution_sequence):
    with screed.open(fasta) as records:
        seq = records["sequence"]
        start = seq.find(evolution_sequence)
        end = start + len(evolution_sequence)
    return start, end


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
@click.option("--insertion-summary", default="insertion_summary.csv", type=click.Path())
@click.option(
    "--insertion-summary-percentages",
    default="insertion_summary_percentages.csv",
    type=click.Path(),
)
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
    insertion_summary,
    insertion_summary_percentages,
    cigar_strings_csv,
    prefix,
):

    if prefix:
        evolved_sequence_counts = f"{prefix}__{evolved_sequence_counts}"
        insertion_summary = f"{prefix}__{insertion_summary}"
        insertion_summary_percentages = f"{prefix}__{insertion_summary_percentages}"
        cigar_strings_csv = f"{prefix}__{cigar_strings_csv}"

    directed_evolution_interval = get_directed_evolution_interval(
        fasta, directed_evolution_sequence
    )

    sequence_counter_df, non_translated_df, cigar_strings_df = count_21nt_7aa(
        bams, directed_evolution_interval
    )
    sequence_counter_df.to_csv(evolved_sequence_counts)
    non_translated_df.to_csv(insertion_summary)
    non_translated_df_percentages = 100 * non_translated_df / non_translated_df.sum()
    non_translated_df_percentages.to_csv(insertion_summary_percentages)

    cigar_strings_df.to_csv(cigar_strings_csv)


if __name__ == "__main__":
    main()
