#!/usr/bin/env python3
""" Generate RNAi insensitive construct """

import argparse
import os
import random
import numpy as np
from typing import NamedTuple, List, TextIO
from Bio.SeqUtils import GC
from Bio.Seq import Seq


# Define our input arguments
class Args(NamedTuple):
    """ Command-line arguments """
    files: List[TextIO]
    out_dir: str
    gc_difference: int
    iterations: int


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Shuffle nucleotides while maintaining amino acids for generating RNAi resistant sequences',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
                        help='Input cDNA',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        nargs='+')

    parser.add_argument('-o',
                        '--out_dir',
                        help='Output directory',
                        metavar='DIR',
                        type=str,
                        default='out')

    parser.add_argument('-g',
                        '--gc_difference',
                        help='Maximum GC content percent difference between input/output',
                        metavar='INT',
                        type=float,
                        nargs='?',
                        default=5)

    parser.add_argument('-i',
                        '--iterations',
                        help='Number of iterations to attempt a good GC match',
                        metavar='INT',
                        type=int,
                        nargs='?',
                        default=1000)

    args = parser.parse_args()

    return Args(args.file, args.out_dir, args.gc_difference, args.iterations)


# --------------------------------------------------
def main() -> None:
    """ Shuffle codons to generate RNAi insensitive constructs """

    args = get_args()

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    num_files, num_seqs = 0, 0
    for fh in args.files:
        num_files += 1
        out_file = os.path.join(args.out_dir, os.path.basename(fh.name))
        out_fh = open(out_file, 'wt')

        # Isolate each sequence in the input file
        for seq in fh:
            num_seqs += 1
            seq = seq.rstrip().upper()

            # Label sequences
            out_fh.write(f'Sequence number: {num_seqs}\n \n')

            # Initial sequence information
            out_fh.write(f'Initial sequence: {seq}\n')
            out_fh.write(f'Starting GC content: {str(round(GC(seq), 2))}\n')
            original_translation = Seq(seq).translate()
            out_fh.write(f'Original translation: {original_translation}\n \n')

            # Isolate codons
            codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]

            # Generate codon optimized
            codon_optimized_list = []

            for codon in codons:
                codon_optimized_list.append(codon_optimize(codon))

            codon_optimized = ''.join(codon_optimized_list)
            codon_optimized_translation = Seq(codon_optimized).translate()
            difference = abs(GC(codon_optimized) - GC(seq))

            out_fh.write(f'Codon optimized: {codon_optimized}\n')
            out_fh.write(
                f'Codon optimized GC content: {str(round(GC(codon_optimized), 2))} (Delta of: {round(difference, 2)})\n')
            out_fh.write(
                f'Codon optimized translation check: {codon_optimized_translation}\n \n')

            # Generate population of candidates with randomized codon usage
            random_codons = []

            if difference > args.gc_difference:
                for i in range(0, args.iterations):
                    temp_random_codons = []
                    for codon in codons:
                        temp_random_codons.append(codon_random(codon))
                    random_codons.append(temp_random_codons)

                # Find one that best matches the GC content of the input string
                random_codon_gc_difference = []
                for randomized_sequence in random_codons:
                    optimal_gc_random_codon_temp = ''.join(randomized_sequence)
                    difference = abs(GC(optimal_gc_random_codon_temp) - GC(seq))
                    random_codon_gc_difference.append(difference)

                minimal_gc_value = min(random_codon_gc_difference)
                minimal_gc_index = random_codon_gc_difference.index(minimal_gc_value)
                optimal_gc_random_codon = ''.join(random_codons[minimal_gc_index])
                codon_randomized_translation = Seq(optimal_gc_random_codon).translate()
                randomized_difference = abs(GC(optimal_gc_random_codon) - GC(seq))

                out_fh.write(f'Closest GC match from {args.iterations} randomized sequences: {optimal_gc_random_codon}\n')
                out_fh.write(
                    f'Codon randomized GC content: {str(round(GC(optimal_gc_random_codon), 2))} (Delta of: {round(randomized_difference, 2)})\n')
                out_fh.write(
                    f'Codon randomized translation check: {codon_randomized_translation}\n \n')

                # Generate per-codon optimization of GC-content
                codon_gc_optimized_list = []

                for codon in codons:
                    codon_gc_optimized_list.append(gc_optimize(codon))

                codon_gc_optimized = ''.join(codon_gc_optimized_list)
                codon_gc_optimized_translation = Seq(codon_gc_optimized).translate()
                gc_optimized_difference = abs(GC(codon_gc_optimized) - GC(seq))

                out_fh.write(f'Per-codon GC optimization sequence: {codon_gc_optimized}\n')
                out_fh.write(
                    f'Per-codon GC optimization GC content: {str(round(GC(codon_gc_optimized), 2))} (Delta of: {round(gc_optimized_difference, 2)})\n')
                out_fh.write(
                    f'Per-codon GC optimization translation check: {codon_gc_optimized_translation}\n')
                out_fh.write('\n -------------------------------------------- \n')
            else:
                out_fh.write('Codon optimized sequence satisfies GC delta requirement')

        out_fh.close()

    print(f'Done, wrote {num_seqs} sequence{"" if num_seqs == 1 else "s"} '
          f'in {num_files} file{"" if num_files == 1 else "s"} '
          f'to directory "{args.out_dir}".')


# --------------------------------------------------
# Codon usage information - source: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=7227

F = [['TTT', .38], ['TTC', .62]]
L = [['TTA', .05], ['TTG', .18], ['CTT', .10], ['CTC', .15], ['CTA', .09], ['CTG', .43]]
I = [['ATT', .34], ['ATC', .47], ['ATA', .19]]
M = [['ATG', 1]]
V = [['GTT', .19], ['GTC', .24], ['GTA', .11], ['GTG', .47]]
S = [['TCT', .08], ['TCC', .24], ['TCA', .09], ['TCG', .20], ['AGT', .14], ['AGC', .25]]
P = [['CCT', .13], ['CCC', .33], ['CCA', .25], ['CCG', .29]]
T = [['ACT', .17], ['ACC', .38], ['ACA', .20], ['ACG', .26]]
A = [['GCT', .19], ['GCC', .45], ['GCA', .17], ['GCG', .19]]
Y = [['TAT', .37], ['TAC', .63]]
stop = [['TAA', .41], ['TAG', .33], ['TGA', .25]]
H = [['CAT', .40], ['CAC', .60]]
Q = [['CAA', .30], ['CAG', .70]]
N = [['AAT', .44], ['AAC', .56]]
K = [['AAA', .30], ['AAG', .70]]
D = [['GAT', .53], ['GAC', .47]]
E = [['GAA', .33], ['GAG', .67]]
C = [['TGT', .29], ['TGC', .71]]
W = [['TGG', 1]]
R = [['CGT', .16], ['CGC', .33], ['CGA', .15], ['CGG', .15], ['AGA', .09], ['AGG', .11]]
G = [['GGT', .21], ['GGC', .43], ['GGA', .29], ['GGG', .07]]

amino_acids = [F, L, I, M, V, S, P, T, A, Y, stop, H, Q, N, K, D, E, C, W, R, G]


# --------------------------------------------------
def codon_optimize(input_sequence: str):
    """" Replace with the most frequently used codon """

    codon_list = []
    frequency_list = []

    aa_opt_key = find_amino_acid_key(input_sequence)

    # Return original sequence if there is no other option
    if len(amino_acids[aa_opt_key]) == 1:
        return input_sequence

    # Make a list of the codons and their relative frequency
    else:
        for codon, frequency in amino_acids[aa_opt_key]:
            codon_list.append(codon)
            frequency_list.append(frequency)

    # Determine which codon is most frequently used
    max_index = frequency_list.index(max(frequency_list))

    # Check if the current value is already the most frequent - if it is remove that value and find second highest
    if codon_list[max_index] == input_sequence:
        codon_list.pop(max_index)
        frequency_list.pop(max_index)
    else:
        pass
    max_index = frequency_list.index(max(frequency_list))

    return codon_list[max_index]


# --------------------------------------------------
def codon_random(input_sequence: str):
    """ Generate a new random codon """

    aa_rand_key = find_amino_acid_key(input_sequence)

    # Return original sequence if there is no other option
    if len(amino_acids[aa_rand_key]) == 1:
        return input_sequence

    # Return a random codon value
    else:
        random_index = random.randint(0, len(amino_acids[aa_rand_key]) - 1)

        while input_sequence == amino_acids[aa_rand_key][random_index][0]:
            random_index = random.randint(0, len(amino_acids[aa_rand_key]) - 1)
        else:
            return amino_acids[aa_rand_key][random_index][0]


# --------------------------------------------------
def gc_optimize(input_sequence: str):
    """ Replaces input codon with the one with the most similar GC content """

    codon_list = []
    gc_list = []

    gc_opt_key = find_amino_acid_key(input_sequence)

    # Return original sequence if there is no other option
    if len(amino_acids[gc_opt_key]) == 1:
        return input_sequence

    # Make a list of the codons and their gc content
    else:
        for codon, _ in amino_acids[gc_opt_key]:
            codon_list.append(codon)
            gc_list.append(GC(codon))

    # Determine which codon has the most similar GC content
    closest_gc_value = closest_value(gc_list, GC(input_sequence))

    # find index of max value
    closest_gc_index = gc_list.index(closest_gc_value)

    # Check if the current value is already the most frequent - if it is remove that value and find second highest
    if codon_list[closest_gc_index] == input_sequence:
        codon_list.pop(closest_gc_index)
        gc_list.pop(closest_gc_index)
    else:
        pass

    closest_gc_value = closest_value(gc_list, GC(input_sequence))
    closest_gc_index = gc_list.index(closest_gc_value)

    return codon_list[closest_gc_index]


# --------------------------------------------------
def closest_value(input_list, input_value):
    """ Find the value in a list closest to the input_value"""
    array = np.asarray(input_list)
    output_value = (np.abs(array - input_value)).argmin()

    return array[output_value]


# --------------------------------------------------
def find_amino_acid_key(input_codon):
    """ Find what amino acid the codon is """
    for index, aa in enumerate(amino_acids):
        for c_index, possible_codons in enumerate(aa):
            if possible_codons[0] == input_codon:
                amino_acid_key = index
                return amino_acid_key


# --------------------------------------------------
if __name__ == '__main__':
    main()
