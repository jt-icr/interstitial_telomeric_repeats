#!/usr/bin/env python3.5

'''Extracts the UCSC genome coordinates from a bed file in a format that
can be entered into the UCSC Genome Browser .'''

from os import listdir


def comma(num):
    return '{:,}'.format(num)  # Takes an int


def process_bed(bedfilename, outfilename, motif_len=30):
    with open(bedfilename, 'r') as fi, open(outfilename, 'w') as fo:
        for line in fi:
            line = line.split()
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            if end - start >= motif_len:
                fo.write('{}:{}-{}\n'.format(chrom, comma(start), comma(end)))


if __name__ == '__main__':

    files = [file for file in listdir('.') if
             file.startswith('uniq') and
             file.endswith('bed')]

    if not files:
        print("No files with extension", file_extension, "found!")

    for file in files:
        outfilename = file + ".coords"
        process_bed(file, outfilename)
