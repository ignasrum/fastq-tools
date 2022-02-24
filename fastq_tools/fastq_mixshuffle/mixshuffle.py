from argparse import ArgumentParser
import os
import traceback
import random
import math

import pysam


def mixshuffle(input_files, output_files, abundance, input_files2=None, seed=None, k=None):
    # check if out filename exists, exit if it does
    out_exists = os.path.isfile(output_files[0])
    if out_exists:
        print(f'Output file already exists: {output_files[0]}')
        return 1
    # check if input is paired
    if input_files2:
        # check if input and input2 are the same length
        if len(input_files) != len(input_files2):
            print(f'Input is paired, but number of files is different: {input_files}, {input_files2}')
            return 1
        # check if two output files are provided
        if len(output_files) != 2:
            print(f'Input is paired, but output files != 2: {output_files}')
            return 1
        out_exists = os.path.isfile(output_files[1])
        if out_exists:
            print(f'Output file already exists: {output_files[1]}')
            return 1
    # check that abundance sums to 100
    if sum(abundance) != 100:
        print(f'Abundance does not sum to 100: {abundance}')
        return 1
    if len(abundance) != len(input_files):
        print(f'Not every input file has an abundance specified: {input_files} {abundance}')
        return 1


    print(f'Input files 1: {input_files}')
    print(f'Input files 2: {input_files2}')
    print(f'Abundance: {abundance}')
    print(f'Output files: {output_files}')
    print(f'Seed: {seed}')
    print(f'k: {k}')

    entries = []
    entries2 = []
    seq_counts = []
    seq_counts2 = []
    choices = []
    try:
        if seed is not None: random.seed(seed)
        # count files/reads
        for in_f, f_idx in zip(input_files, range(len(input_files))):
            entries.append([])
            choices.append(f_idx)
            with pysam.FastxFile(in_f) as f:
                i = 0
                for entry in f:
                    entries[f_idx].append(entry)
                    i += 1
                seq_counts.append(i)
        if input_files2:
            for in_f, f_idx in zip(input_files2, range(len(input_files2))):
                entries2.append([])
                with pysam.FastxFile(in_f) as f:
                    i = 0
                    for entry in f:
                        entries2[f_idx].append(entry)
                        i += 1
                    seq_counts2.append(i)
            # check if number of reads is the same in input files 
            for c1, c2 in zip(seq_counts, seq_counts2):
                if c1 != c2:
                    print(f'Input files are paired, but number of reads is different: {seq_counts} != {seq_counts2}')
                    return 1

        if k == None: k = sum(seq_counts)

        if k is not None and k > sum(seq_counts):
            print(f'k is larger than total reads: {k} > {seq_counts}')
            return 1

        # check if enough reads are present in input files to satisfy abundance requirements
        if seq_counts[0] < k*(abundance[0]/100):
            print(f'Not enough reads in input file 1 to satisfy abundance requirements.')
            print(f'Reads in input file 1: {seq_counts[0]}')
            print(f'Required reads: {k*(abundance[0]/100)}')
            return 1
        if seq_counts[1] < k*(abundance[1]/100):
            print(f'Not enough reads in input file 2 to satisfy abundance requirements.')
            print(f'Reads in input file 2: {seq_counts[1]}')
            print(f'Required reads: {k*(abundance[1]/100)}')
            return 1

        # calculate which sequences from each file will be picked
        final_choices = []
        for f in range(len(input_files)):
            pop = [(f, x) for x in range(0, seq_counts[f]-1)]
            new_k = math.floor(k * (abundance[f] / 100))
            print(f'len(pop): {len(pop)}')
            print(f'new_k: {new_k}')
            final_choices += random.sample(pop, new_k)
        random.shuffle(final_choices)

        # write these sequences to output files
        with open(output_files[0], mode='w') as fout:
            for choice in final_choices:
                entry = entries[choice[0]][choice[1]]
                # write sequence name
                fout.write(f'@{entry.name}\n')
                # write sequence
                fout.write(f'{entry.sequence}\n')
                # write quality line break
                fout.write(f'+\n')
                # write nucleotide quality
                fout.write(f'{entry.quality}\n')
        if input_files2:
            # write these sequences to output files
            with open(output_files[1], mode='w') as fout:
                for choice in final_choices:
                    entry = entries2[choice[0]][choice[1]]
                    # write sequence name
                    fout.write(f'@{entry.name}\n')
                    # write sequence
                    fout.write(f'{entry.sequence}\n')
                    # write quality line break
                    fout.write(f'+\n')
                    # write nucleotide quality
                    fout.write(f'{entry.quality}\n')
    except Exception as e:
        print(f'Exception encoutered while reading/writing FastQ files: {e}')
        print(traceback.format_exc())
        return 1
    return 0

def main():
    parser = ArgumentParser(description='Mix and/or shuffle FastQ files')

    parser.add_argument('-1',
                        '--input',
                        type=str,
                        required=True,
                        nargs=('*'),
                        metavar=('<in.fastq>'),
                        help='str: path to input FastQ files, count has to be the same length as --input2')
    parser.add_argument('-2',
                        '--input2',
                        type=str,
                        nargs=('*'),
                        metavar=('<in.fastq>'),
                        help='str: path to paired input FastQ files, count has to be the same length as --input')
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        required=True,
                        nargs=('*'),
                        metavar=('<out.fastq>'),
                        help='str: path to output FastQ file (count has to be two if --input2 is provided)')
    parser.add_argument('-a',
                        '--abundance',
                        type=int,
                        nargs=('*'),
                        metavar=('<int>'),
                        help='int: percentage of reads from each FastQ file, has to sum to 100, n numbers if n input files')
    parser.add_argument('-s',
                        '--seed',
                        type=int,
                        metavar=('<int>'),
                        help='int: seed for random number generator, default: random')
    parser.add_argument('-k',
                        type=int,
                        metavar=('<int>'),
                        help='int: number of reads in output FastQ, cannot be more than total of both FastQ files, default: sum of reads in all input files')

    args = parser.parse_args()
    input_files = args.input
    input_files2 = args.input2
    output_files = args.output
    abundance = args.abundance
    seed = args.seed
    k = args.k

    result = mixshuffle(input_files, output_files, abundance, input_files2=input_files2, seed=seed, k=k)
    print(f'Exited with result: {result}')
    return result

if __name__ == '__main__':
    main()
