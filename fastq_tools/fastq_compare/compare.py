import itertools
from argparse import ArgumentParser

import pysam


def count(file, separator):
    entries = {}
    with pysam.FastxFile(file) as f:
        for record in f:
            seq = record.name.split(separator)[0]
            if seq not in entries:
                entries[seq] = 1
            else: entries[seq] += 1
    total = 0
    for entry in entries:
        total += entries[entry]
    return entries, total


def print_cmp(msg, file1_num, file2_num):
    print(msg)
    print('\tfile1: ', file1_num)
    print('\tfile2: ', file2_num)
    print('\tdifference: ', file1_num - file2_num)
    print('\tdifference (%): ', ((file1_num - file2_num) / file1_num) * 100)

def compare(file1, file2, separator):
    entries1, total1 = count(file1, separator)
    entries2, total2 = count(file2, separator)
    for entry1, entry2 in itertools.zip_longest(entries1, entries2):
        if entry1 is not None:
            val = entries2[entry1] if entry1 in entries2 else 0
            print_cmp(entry1, entries1[entry1], val)
        if entry2 is not None and entry2 not in entries1:
            print_cmp(entry2, 0, entries2[entry2])
    print('---------------')
    print_cmp('total reads: ', total1, total2)

def main():
    parser = ArgumentParser(description='Compare two FastQ files')

    parser.add_argument('file1', type=str,
                        help='str: path to file 1')
    parser.add_argument('file2', type=str,
                        help='str: path to file 2')
    parser.add_argument('separator', type=str,
                        help='char: separator used when splitting sequence names')

    args = parser.parse_args()
    file1 = args.file1
    file2 = args.file2
    separator = args.separator

    compare(file1, file2, separator)

if __name__ == '__main__':
    main()
