import argparse
import os
import re

__author__ = 'oyasnev'

QUAST_CONTIGS = "/contigs_reports/contigs_report_contigs.stdout"

FORWARD_STRAND = "FS"
REVERSE_COMPLEMENT_STRAND = "RC"


class RealAlignment:
    def __init__(self, align_tuple):
        self.ref_pos1 = int(align_tuple[0])
        self.ref_pos2 = int(align_tuple[1])
        self.contig_pos1 = int(align_tuple[2])
        self.contig_pos2 = int(align_tuple[3])
        if self.contig_pos2 - self.contig_pos1 >= 0:
            self.strand = FORWARD_STRAND
        else:
            self.strand = REVERSE_COMPLEMENT_STRAND


class ExtMisassembly:
    def __init__(self):
        self.contig_name = ''
        self.type = ''
        self.align_list = []


""" Parse contigs with extensive misassemblies """


def parse_contigs(contigs_path):
    print("Parsing contigs with extensive misassemblies...")
    try:
        fp = open(contigs_path, "r")
    except:
        print("ERROR: some error on reading {}".format(contigs_path))
        exit(1)
    else:
        with fp:
            # At first, we find paragraphs with an extensive misassembly
            # Pattern:
            #   "CONTIG" + several lines + a line with "Extensive misassembly" +
            #   + several lines + empty line.
            text = fp.read()
            contig_pattern = re.compile("CONTIG(?:.+\n)+.+Extensive misassembly.+\n(?:.+\n)+\n")

            # Patterns for misassembly info:
            # Contig name: "NODE_" + some characters until a whitespace
            node_pattern = re.compile("NODE_[^\s]+")
            # Misassembly info block: "Real Alignment" + some symbols and EOL +
            #   + "Extensive misassembly" + symbols and EOL + another "Real Alignment"
            # (the pattern supports overlapping)
            mis_pattern = re.compile("(?=(Real Alignment .+\n.+Extensive misassembly.+\n.+Real Alignment.+\n))")
            # Real alignment: "Real alignment %number%: %ref_pos1% %ref_pos2% | %contig_pos1% %contig_pos2%"
            align_pattern = re.compile("Real Alignment \d+: (\d+) (\d+) \| (\d+) (\d+)")
            # Misassembly type: "Extensive misassembly ( %type%"
            mis_type_pattern = re.compile("Extensive misassembly \( (\w+)")

            # Parse contig blocks
            mis_list = []
            for contig_text in contig_pattern.findall(text):
                contig_name = node_pattern.search(contig_text).group(0)
                for mis_text in mis_pattern.findall(contig_text):
                    mis = ExtMisassembly()
                    mis.contig_name = contig_name
                    # Real alignment
                    for align_tuple in align_pattern.findall(mis_text):
                        mis.align_list.append(RealAlignment(align_tuple))
                    # Misassembly type
                    mis.type = mis_type_pattern.search(mis_text).group(1)
                    mis_list.append(mis)
    print("Done")
    print("{} extensive misassemblies found".format(len(mis_list)))
    print()
    return mis_list


""" Parse and validate input arguments """


def get_args():
    parser = argparse.ArgumentParser(description='Miss Classification.')
    parser.add_argument("quast", help="path to QUAST report folder")
    args = parser.parse_args()

    print("Validating input arguments...")
    args.quast_contigs = args.quast + QUAST_CONTIGS
    if not os.path.exists(args.quast_contigs):
        print("ERROR: {} is not a correct QUAST report folder".format(args.quast))
        exit(1)

    def test_file(path):
        try:
            fp = open(path, "r")
        except IOError as e:
            print("ERROR: {} on reading {}".format(e.strerror, args.quast_contigs))
            exit(1)
        except:
            print("ERROR: unknown error on reading {}".format(args.quast_contigs))
            exit(1)
        else:
            fp.close()

    test_file(args.quast_contigs)

    print("OK")
    print()
    return args


#############################
# MAIN

args = get_args()
mis_list = parse_contigs(args.quast_contigs)
