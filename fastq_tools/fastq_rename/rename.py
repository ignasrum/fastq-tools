import itertools
from argparse import ArgumentParser
import os
import traceback

import pysam


def rename(input_file, output_file, base_name, append_end=None):
    append_end = append_end if append_end != None else ''
    try:
        i = 1
        with pysam.FastxFile(input_file) as fin, open(output_file, mode='w') as fout:
            for entry in fin:
                # write sequence name
                fout.write(f'@{base_name}{i}{append_end}\n')
                # write sequence
                fout.write(f'{entry.sequence}\n')
                # write quality line break
                fout.write(f'+\n')
                # write nucleotide quality
                fout.write(f'{entry.quality}\n')
                i += 1
    except Exception as e:
        print(f'Exception encoutered while reading/writing FastQ files: {e}')
        print(traceback.format_exc())
        return 1
    return 0

def main():
    parser = ArgumentParser(description='Rename FastQ files')

    parser.add_argument('-i',
                        '--input',
                        type=str,
                        required=True,
                        metavar=('<in.fastq>'),
                        help='str: path to input FastQ file')
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        required=True,
                        metavar=('<out.fastq>'),
                        help='str: path to output FastQ file')
    parser.add_argument('-b',
                        '--base',
                        type=str,
                        required=True,
                        metavar=('<str>'),
                        help='str: new base name for sequences in input file')
    parser.add_argument('-e',
                        '--end',
                        type=str,
                        metavar=('<str>'),
                        help='str: append this to end of sequence names in input file')

    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    base_name = args.base
    append_end = args.end

    # check if out filename exists, exit if it does
    out_exists = os.path.isfile(output_file)
    if out_exists:
        print(f'Output file already exists: {output_file}')
        return 1

    result = rename(input_file, output_file, base_name, append_end)
    print(f'Exited with result: {result}')

if __name__ == '__main__':
    main()
