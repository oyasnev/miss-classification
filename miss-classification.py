import argparse
from enum import Enum
import os
import re

__author__ = 'oyasnev'

QUAST_CONTIGS = "/contigs_reports/contigs_report_contigs.stdout"

DISTANCE_ESTIMATION_INPUT_FILE = "/__miss_classification_dist_est_input.txt"

OVERLAP_THRESHOLD = 200


class Strand(Enum):
    forward = "fs"
    reverseComplement = "rc"


class RealAlignment:
    def __init__(self, align_tuple):
        self.ref_pos1 = int(align_tuple[0])
        self.ref_pos2 = int(align_tuple[1])
        self.contig_pos1 = int(align_tuple[2])
        self.contig_pos2 = int(align_tuple[3])
        if self.contig_pos2 - self.contig_pos1 >= 0:
            self.strand = Strand.forward
        else:
            self.strand = Strand.reverseComplement


class ExtMisassembly:
    def __init__(self):
        self.contig_name = ''
        self.type = ''
        self.align_list = []

    def overlap_length(self):
        pos1 = max(self.align_list[0].contig_pos1, self.align_list[0].contig_pos2)
        pos2 = min(self.align_list[1].contig_pos1, self.align_list[1].contig_pos2)
        if pos2 > pos1:
            return 0
        else:
            return pos1 - pos2


class MisClassification:
    def __init__(self):
        self.broken_bone = []
        self.ignored = []
        self.unknown = []


""" Predict misassemblies classification """


def predict_classes(mis_list):
    mis_class = MisClassification()
    for mis in mis_list:
        # broken bone
        if mis.align_list[0].strand is Strand.forward and mis.align_list[1].strand is Strand.forward:
            if mis.overlap_length() >= OVERLAP_THRESHOLD:
                mis_class.broken_bone.append(mis)
            else:
                mis_class.ignored.append(mis)
            continue

        # other cases
        mis_class.unknown.append(mis)

    print("Predicted classification")
    print("Broken bone: {}".format(len(mis_class.broken_bone)))
    print("Ignored:     {}".format(len(mis_class.ignored)))
    print("Unknown:     {}".format(len(mis_class.unknown)))
    print()

    return mis_class


""" Write info for distance estimation stage """


def write_for_dist_est(mis_class, args):
    print("Writing info for distance estimation stage...")
    try:
        fp = open(args.spades + DISTANCE_ESTIMATION_INPUT_FILE, "w")
    except:
        print("ERROR: some error on writing {}".format(args.spades + DISTANCE_ESTIMATION_INPUT_FILE))
        exit(1)
    else:
        with fp:
            fp.write(args.contigs + "\n")
            fp.write(args.ref + "\n")
            # broken bone
            fp.write(str(len(mis_class.broken_bone)))
            for mis in mis_class.broken_bone:
                fp.write("{} {} {} {}\n".format(
                    mis.align_list[0].ref_pos1, mis.align_list[0].ref_pos2,
                    mis.align_list[0].contig_pos1, mis.align_list[0].contig_pos2))
                fp.write("{} {} {} {}\n".format(
                    mis.align_list[1].ref_pos1, mis.align_list[1].ref_pos2,
                    mis.align_list[1].contig_pos1, mis.align_list[1].contig_pos2))


            print("Done")
            print()


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
    parser.add_argument("spades", help="path to SPAdes folder")
    parser.add_argument("contigs", help="path to contigs FASTA file")
    parser.add_argument("ref", help="path to reference FASTA file")
    args = parser.parse_args()

    print("Validating input arguments...")

    def is_file_readable(path):
        if not os.path.exists(path):
            print("ERROR: {} does not exits".format(path))
            exit(1)
        try:
            fp = open(path, "r")
        except IOError as e:
            print("ERROR: {} on reading {}".format(e.strerror, path))
            exit(1)
        except:
            print("ERROR: unknown error on reading {}".format(path))
            exit(1)
        else:
            fp.close()

    def is_folder_writable(path):
        try:
            test_file_path = path + "/__miss_classification_test_file.txt"
            fp = open(test_file_path, "w")
        except IOError as e:
            print("ERROR: {} on writing in {}".format(e.strerror, path))
            exit(1)
        except:
            print("ERROR: unknown error on writing {}".format(path))
        else:
            fp.close()
            os.remove(test_file_path)

    args.quast_contigs = args.quast + QUAST_CONTIGS
    if not os.path.exists(args.quast_contigs):
        print("ERROR: {} is not a correct QUAST report folder".format(args.quast))
        exit(1)
    is_file_readable(args.quast_contigs)

    if not os.path.exists(args.spades):
        print("ERROR: {} is not a correct SPAdes folder".format(args.spades))
        exit(1)
    is_folder_writable(args.spades)

    is_file_readable(args.contigs)
    is_file_readable(args.ref)

    print("OK")
    print()
    return args


#############################
# MAIN

args = get_args()
mis_list = parse_contigs(args.quast_contigs)

mis_class = predict_classes(mis_list)

write_for_dist_est(mis_class, args)
