#!/usr/bin/env python3.5

'''This module was developed under python3.5. It contains two classes
and several module level methods for the following two purposes.

1) Find both forward and reverse telomere repeats, including their most
common degenerate forms (TT.GGG and CCC.AA). It can also make a bed file
with the telomere genome coordinates for a chromosome that can be uploaded
in a genome browser or used for additional analyses. This part of the module
can be run as a batch script using the batch_fasta_telo() method in a
directory of fasta chromosome files. Output will print to stdout by default
or can be directed to file, or both with tee. Bed files will be automatically
generated in the same directory for each chromosome.

2) Find the intersects of two bed files, primarily the intersect of a
telomere bed file with another bed file of some ENCODE feature type. Also
outputs the intersect of both to bedfile format, allows for merging, removal
of duplicates, and counting of intersects.

All code in this module conforms to the PEP 0008 Python style guide with
the aid of the pep8 program.

Dependancy: pybedtools [https://daler.github.io/pybedtools/]'''

from os import listdir, chdir
from shutil import copyfileobj
import io
import re
import time

from pybedtools import BedTool

__author__ = "Jeffrey P Tomkins, PhD"
__contributor__ = ""
__copyright__ = "Copyright 2015, Institute for Creation Research"
__credits__ = []
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey P Tomkins"
__email__ = "jtomkins@icr.org"
__status__ = "Production"

rm_nondna = re.compile(' |>.*$|\d|\t|\n')


class TeloFinder:

    def __init__(self, fasta_filename):
        '''Creates an instance of TeloFinder from a fasta format file which
        contains the results of a search of both forward and reverse telomere
        repeats as described in the module documentation. Also generates their
        index positions corresponding to genome coordinates using the UCSC
        genome coordinate index format. NOTE! Index will be off by one for
        ensembl!'''

        def _from_fasta(fasta_infilename):
            '''Get the DNA seq from a FASTA file as one contiguous string.'''
            try:
                with io.open(fasta_infilename, 'r', encoding='utf8') as fi:
                    return ''.join([rm_nondna.sub('', line.upper())
                                    for line in fi])
            except IOError as e:
                print('Cannot open', dna_file, e.errno, e.strerror)

        dna = _from_fasta(fasta_filename)

        self.exact_fhits = [(m.start(0), m.end(0))
                            for m in re.finditer(r'((TTAGGG){2,})', dna)]
        self.exact_rhits = [(m.start(0), m.end(0))
                            for m in re.finditer(r'((CCCTAA){2,})', dna)]
        self.inexact_fhits = [(m.start(0), m.end(0))
                              for m in re.finditer(r'((TT.GGG){2,})', dna)]
        self.inexact_rhits = [(m.start(0), m.end(0))
                              for m in re.finditer(r'((CCC.AA){2,})', dna)]

        # Results from re.finditer above
        self.telo_results = (self.exact_fhits, self.exact_rhits,
                             self.inexact_fhits, self.inexact_rhits)

    def summary(self):
        '''Print the results with descriptive headers'''

        def _mer_num(beg, end):
            '''This Function finds the length of each hit and divides by
            6 bases for number of telos in repeat.'''
            return (end - beg) // 6

        def _hit_summary(hits):
            '''Hits are returned from re.findall in a list of 2-element tuples.
            This fucntion uses a dict to tally the number of hits per telo
            repeat. Returns a list of the tallies as strings.'''
            hit_list = [_mer_num(hit[0], hit[1]) for hit in hits]
            hit_tallies = {x: hit_list.count(x) for x in hit_list}
            return ['{}_mer => {}'.format(mer, hit_num)
                    for mer, hit_num in hit_tallies.items()]

        telo_varnames = ['exact_fhits', 'exact_rhits',
                         'inexact_fhits', 'inexact_rhits']

        for telo_varname, telo_result in zip(telo_varnames, self.telo_results):
            print(telo_varname)
            for item in _hit_summary(telo_result):
                print(item)

    def make_bed(self, chr):
        '''Creates a bed file from re.finditer data in format below.
        chr22   11707814    11707825    forward_telo    0   +
        '''
        # Names for column 4 in bed file
        names = ['forward_telo', 'reverse_telo',
                 'degen_forward_telo', 'degen_reverse_telo']

        try:
            with open("telos_{}.bed".format(chr), 'w') as fo:
                # Combine names and results in a dict and make bed file
                telo_dict = dict(zip(names, self.telo_results))
                for name, result in telo_dict.items():
                    for item in result:
                        fo.write("{}\t{}\t{}\t{}\t0\t+\n".format(
                                chr.replace('_', ''), item[0], item[1], name))
        except IOError as e:
            print('Cannot open file to write', e.errno, e.strerror)


class TeloIntersect:

    def __init__(self, telo_bedfile, other_bedfile):
        '''Create an instance of intersections between the telo repeats and
        some other feature in bed file format. Methods provide a way to
        print to file the intersected bed file. Code brevity made possible
        by the pybedtools module (Dale, Pederson and Quinlan (2011),
        Bioinformatics. 27:3423-3424). More functions to be added to this
        class in the future.'''

        self.telo_bed = BedTool(telo_bedfile)
        self.other_bed = BedTool(other_bedfile)

    def intersect_bed(self, outbed_name='intersect_out.bed'):
        self.telo_and_other = self.telo_bed.intersect(self.other_bed)
        self.telo_and_other.saveas(outbed_name)


def batch_fasta_telo():
    '''Enables this module to be run as a batch script - processing all
    the FASTA files in the directory in which the script is located.'''

    # Put chr files in list - formats: *.01*.fa,...*.22*.fa,*._X*.fa, etc
    chr_files = sorted([file for file in listdir('.')
                        if file.endswith('fa')])

    if not chr_files:
        print("No files with *.fa found!")

    # Extract the chrom name from each fasta filename and put in list
    chrs = ['chr' + chr.rsplit('me.', 1)[1][:2] for chr in chr_files]

    # Make an iter and remove the chr '0' and '_' for proper bed files
    chr = iter([chr.replace('r0', 'r').replace('_', '') for chr in chrs])

    # Get the telos, their counts, and print report for each chrom
    for chr_file in chr_files:
        obj = TeloFinder(chr_file)
        print("<{}>".format(chr_file))
        obj.summary()
        print()
        # Make telo bed files for each chrom
        obj.make_bed(next(chr))


def concat_beds(out_bed='all.bed'):
    '''Concatenates all bed files in directory - outfile as arg.'''

    bed_files = sorted([file for file in listdir('.')
                        if file.endswith('bed')])

    if not bed_files:
        print("No files with *.bed found!")

    try:
        fo = open(out_bed, 'w')
    except IOError:
        print('Cannot open', out_bed)

    for bed_file in bed_files:
        try:
            fi = open(bed_file, 'r')
            copyfileobj(fi, fo)
            fi.close()
  bed_list      except IOError as e:
            print('Cannot open', bed_file, e.errno, e.strerror)
    fo.close()


def uniq_bed(bed_file, out_bedfile='uniq.bed'):
    '''Remove duplicates in bedfile and print to new file
    maintaining order of bed entries.'''
    try:
        with open(bed_file, 'r') as fi:
            bed_list = [line for line in fi]
            perfect_list = []
            degen_list = []
            uniq_bed = []
            # Separate perfect and degen telos
            for item in bed_list:
                if item.split()[3].startswith('degen'):
                    degen_list.append(item)
                else:
                    perfect_list.append(item)
            # Eliminate degen telo entries that are also perfect
            tmp_degen_list = [item.split()[0:3] for item in degen_list]
            for p_item in perfect_list:
                if p_item.split()[0:3] in tmp_degen_list:
                    if item.split()[3].startswith('degen') in bed_list:
                        bed_list.remove(item)
            # Remove duplicates
            if item not in uniq_bed:
                    uniq_bed.append(item)
    except IOError as e:
            print('Cannot open', bed_file, e.errno, e.strerror)

    try:
        with open(out_bedfile, 'w') as fo:
            for item in uniq_bed:
                fo.write(item)
    except IOError as e:
            print('Cannot open', out_file, e.errno, e.strerror)


def ct_bed_telos(bed_file):
    try:
        with open(bed_file, 'r') as fi:
            ttaggg = tt_ggg = ccctaa = ccc_aa = linect = 0
            for line in fi:
                ttaggg += line.count('\tforward_telo')
                tt_ggg += line.count('degen_forward_telo')
                ccctaa += line.count('\treverse_telo')
                ccc_aa += line.count('degen_reverse_telo')
                linect += 1
            print("<{}>".format(bed_file))
            print("ttaggg:", ttaggg)
            print("tt_ggg:", tt_ggg)
            print("ccctaa:", ccctaa)
            print("ccc_aa:", ccc_aa)
            print("linect:", linect)
    except IOError as e:
        print('Cannot open', bed_file, e.errno, e.strerror)


if __name__ == '__main__':
    # Get the telo bedfiles and stats for each chrom
    # Pipe stdout with tee to screen and file
    batch_fasta_telo()
    # Put the concatenated output bed file in separate directory with
    # the hg38 feature bed files to be used for intersect analyses
    concat_beds('bed_analysis/allchroms_telo_hg38.bed')
    # Do intersect analyses in separate directory
    chdir('bed_analysis/')
    tfile = 'allchroms_telo_hg38.bed'

    # Get intersect with gencode v22
    telo_gc22 = TeloIntersect(tfile, 'gencode22_hg38.bed')
    telo_gc22.intersect_bed('intersect_telo_gc22.bed')

    # Get intersect with UCSC lincRNAs
    telo_ucsclinc = TeloIntersect(tfile, 'human_lincRNAs_hg38_ucsc.bed')
    telo_ucsclinc.intersect_bed('intersect_telo_ucsclinc.bed')

    # Get intersect with Hangauer lincRNAs
    telo_hangauerlinc = TeloIntersect(tfile, 'hangauer_lincs_hg38.bed')
    telo_hangauerlinc.intersect_bed('intersect_telo_hangauerlinc.bed')

    # Get intersect with ReMap TF binding peaks
    telo_remap = TeloIntersect(tfile,
                               'nrPeaks_TFbinding_ReMap_all_merged_hg38.bed')
    telo_remap.intersect_bed('intersect_telo_remap_tfbinding.bed')

    # Get intersect with FANTOM5 enhancer-tss sites
    telo_enhancer_tss = TeloIntersect(tfile,
                                      'enhancer_tss_associations_hg38.bed')
    telo_enhancer_tss.intersect_bed('intersect_telo_enhancer_tss.bed')

    # Get intersect with FANTOM5 permissive enhancers
    telo_permissive_enhancer = TeloIntersect(tfile,
                                             'permissive_enhancers_hg38.bed')
    telo_permissive_enhancer.intersect_bed(
                                    'intersect_telo_permissive_enhancer.bed')

    # Get intersect with FANTOM5 robust enhancers
    telo_robust_enhancers = TeloIntersect(tfile, 'robust_enhancers_hg38.bed')
    telo_robust_enhancers.intersect_bed('intersect_telo_robust_enhancers.bed')

    # Get intersect with FANTOM5 ubiquitous enhancers
    telo_ubiquitous_enhancers = TeloIntersect(tfile,
                                              'ubiquitous_enhancers_hg38.bed')
    telo_ubiquitous_enhancers.intersect_bed(
                                    'intersect_telo_ubiquitous_enhancers.bed')

    telo_trusight_inherited_disease = TeloIntersect(tfile,
    						'trusight_inherited_disease_manifest_a_hg38.bed')
    telo_trusight_inherited_disease.intersect_bed(
    	                     'intersect_telo_trusight_inherited_disease.bed')

    # Remove duplicate lines in intersected bed files with prefix 'uniq'
    # in new filename
    for file in listdir('.'):
        if file.startswith('intersect'):
            newfile = 'uniq_' + file
            uniq_bed(file, newfile)

    # Print number of entries in each uniq*.bed file
    for file in listdir('.'):
        if file.startswith('uniq'):
            ct_bed_telos(file)
