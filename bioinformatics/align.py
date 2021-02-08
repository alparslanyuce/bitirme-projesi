from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio import Align
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.FastaIO import FastaIterator
from Bio import pairwise2
from Bio.SearchIO import BlatIO
import subprocess
import enum
import timeit


class Algorithm(enum.Enum):
    Blast = 1,
    Blat = 2,
    PSA = 3,
    SOAP2 = 4,
    PatternHunter = 5


class Aligner:
    def __init__(self, db):
        self._db = db
        self._E_VALUE_THRESH = 1e-20
        self._LINE_LENGTH = 100

    def _get_html(self, filename, score, query, match, subject):
        html = "File: {}".format(filename) + '</br>'
        html += "Score: {}".format(score) + '</br>'
        html += '<pre>'

        html += ''.join('-' for i in range(self._LINE_LENGTH)) + '\n'

        for i in range(0, len(query), self._LINE_LENGTH):
            html += query[i:i + self._LINE_LENGTH] + '\n'
            html += match[i:i + self._LINE_LENGTH] + '\n'
            html += subject[i:i + self._LINE_LENGTH] + '\n'

            if i + self._LINE_LENGTH < len(query):
                html += '\n'

        html += ''.join('-' for i in range(self._LINE_LENGTH)) + '\n\n'

        html += '</pre>'

        return html

    def align(self, sequence_file, algorithm):
        if algorithm is Algorithm.Blast:
            return self._blast(sequence_file)
        elif algorithm is Algorithm.Blat:
            return self._blat(sequence_file)
        elif algorithm is Algorithm.PSA:
            return self._psa(sequence_file)

    def _blast(self, sequence_file):
        html = ''
        output_file = 'results.xml'

        db = open(self._db, 'r')
        start = timeit.default_timer()

        for filename in db.readlines():
            blastn_cmd = NcbiblastnCommandline(query=sequence_file, subject=filename, outfmt=5, out=output_file)
            stdout, stderr = blastn_cmd()

            for record in NCBIXML.parse(open(output_file)):
                if record.alignments:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect < self._E_VALUE_THRESH:
                                html += self._get_html(filename, hsp.score, hsp.query, hsp.match, hsp.sbjct)

        stop = timeit.default_timer()
        html += 'Time elapsed: {}'.format(stop - start)

        return html

    def _blat(self, sequence_file):
        html = ''
        output_file = 'results'

        db = open(self._db, 'r')
        start = timeit.default_timer()

        for filename in db.readlines():
            blat_cmd = ['blat', sequence_file, filename.strip(), 'out=blast', output_file]
            subprocess.run(blat_cmd)

            try:
                for record in SearchIO.read(output_file, 'blast-text'):
                    for hsp in record.hsps:
                        html += self._get_html(filename, hsp.aln_span, str(hsp.query.seq).upper(),
		                                       hsp.aln_annotation['similarity'], str(hsp.hit.seq).upper())
            except:
                continue

        stop = timeit.default_timer()
        html += 'Time elapsed: {}'.format(stop - start)

        return html

    def _psa(self, sequence_file): # Performs pairwise sequence alignment using dynamic programming
        # db = open(self._db, 'r')
        # result = ""
        #
        # query_sequence = ""
        #
        # with open(sequence_file) as file_handle:
        #     for iterator in FastaIterator(file_handle):
        #         query_sequence += iterator.seq
        #
        # for file in db.readlines():
        #     with open(file.strip()) as file_handle:
        #         for iterator in FastaIterator(file_handle):
        #             for alignment in pairwise2.align.localxx(query_sequence, iterator.seq):
        #                result += str(format_alignment(*alignment))

        html = ''

        query_sequence = ''

        with open(sequence_file) as file_handle:
            for iterator in FastaIterator(file_handle):
                query_sequence += iterator.seq

        filename = 'db/sequence-psa.fasta'

        start = timeit.default_timer()

        with open(filename) as file_handle:
            for iterator in FastaIterator(file_handle):
                for alignment in pairwise2.align.globalxx(query_sequence, iterator.seq):
                    query, match, subject = format_alignment(*alignment)

                    html += self._get_html(filename, alignment.score, query, match, subject)

        stop = timeit.default_timer()
        html += 'Time elapsed: {}'.format(stop - start)

        return html


def format_alignment(align1, align2, score, begin, end, full_sequences=False):
    align_begin = begin
    align_end = end
    start1 = start2 = ""
    start_m = begin  # Begin of match line (how many spaces to include)
    # For local alignments:
    if not full_sequences and (begin != 0 or end != len(align1)):
        # Calculate the actual start positions in the un-aligned sequences
        # This will only work if the gap symbol is '-' or ['-']!
        start1 = str(len(align1[:begin]) - align1[:begin].count("-") + 1) + " "
        start2 = str(len(align2[:begin]) - align2[:begin].count("-") + 1) + " "
        start_m = max(len(start1), len(start2))
    elif full_sequences:
        start_m = 0
        begin = 0
        end = len(align1)

    if isinstance(align1, list):
        # List elements will be separated by spaces, since they can be
        # of different lengths
        align1 = [a + " " for a in align1]
        align2 = [a + " " for a in align2]

    s1_line = ["{:>{width}}".format(start1, width=start_m)]  # seq1 line
    m_line = [" " * start_m]  # match line
    s2_line = ["{:>{width}}".format(start2, width=start_m)]  # seq2 line

    for n, (a, b) in enumerate(zip(align1[begin:end], align2[begin:end])):
        # Since list elements can be of different length, we center them,
        # using the maximum length of the two compared elements as width
        m_len = max(len(a), len(b))
        s1_line.append("{:^{width}}".format(a, width=m_len))
        s2_line.append("{:^{width}}".format(b, width=m_len))
        if full_sequences and (n < align_begin or n >= align_end):
            m_line.append("{:^{width}}".format(" ", width=m_len))  # space
            continue
        if a == b:
            m_line.append("{:^{width}}".format("|", width=m_len))  # match
        elif a.strip() == "-" or b.strip() == "-":
            m_line.append("{:^{width}}".format(" ", width=m_len))  # gap
        else:
            m_line.append("{:^{width}}".format(".", width=m_len))  # mismatch

    return "".join(s1_line)[6:], "".join(m_line)[6:], "".join(s2_line)[6:]
