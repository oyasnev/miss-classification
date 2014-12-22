import argparse
from enum import Enum
import os
import re
import subprocess

__author__ = 'oyasnev'

QUAST_CONTIGS = "/contigs_reports/contigs_report_contigs.stdout"

DISTANCE_ESTIMATION_INPUT_FILENAME = "/__miss_classification_dist_est_input.txt"
DISTANCE_ESTIMATION_OUTPUT_FILENAME = "__miss_classification_dist_est_output.txt"
SIMPLIFICATION_OUTPUT_FILENAME = "__miss_classification_simpl_output.txt"


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

    def __str__(self):
        return "{}\n  ref1: [{} - {}] --> [{} - {}]\n  ref2: [{} - {}] --> [{} - {}]".format(
            self.contig_name,
            self.align_list[0].ref_pos1, self.align_list[0].ref_pos2, self.align_list[0].contig_pos1, self.align_list[0].contig_pos2,
            self.align_list[1].ref_pos1, self.align_list[1].ref_pos2, self.align_list[1].contig_pos1, self.align_list[1].contig_pos2
        )

    def overlap_length(self):
        pos1 = max(self.align_list[0].contig_pos1, self.align_list[0].contig_pos2)
        pos2 = min(self.align_list[1].contig_pos1, self.align_list[1].contig_pos2)
        if pos2 > pos1:
            return 0
        else:
            return pos1 - pos2


class MisClassification:
    def __init__(self):
        # What we predict based on Quast
        self.broken_bone = []
        self.ignored = []
        self.unknown = []
        # What Spades says
        self.spades_broken_bone = []


class BrokenBoneSimplOutput:
    def __init__(self):
        self.is_ready = False
        self.is_start = False
        self.is_start_broken = False
        self.start_coverage = 0.0
        self.is_end = False
        self.is_end_broken = False
        self.end_coverage = 0.0


def predict_classes(mis_list):
    """ Predict misassemblies classification """
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


def write_dist_est_input(mis_class, args):
    """ Write input data for distance estimation stage """
    print("Writing input data for distance estimation stage...")
    try:
        fp = open(args.spades + DISTANCE_ESTIMATION_INPUT_FILENAME, "w")
    except:
        print("ERROR: some error on writing {}".format(args.spades + DISTANCE_ESTIMATION_INPUT_FILENAME))
        exit(1)
    else:
        with fp:
            fp.write(args.contigs + "\n")
            fp.write(args.ref + "\n")
            # broken bone
            fp.write(str(len(mis_class.broken_bone)) + "\n")
            for mis in mis_class.broken_bone:
                fp.write(mis.contig_name + "\n")
                fp.write("{} {} {} {}\n".format(
                    mis.align_list[0].ref_pos1, mis.align_list[0].ref_pos2,
                    mis.align_list[0].contig_pos1, mis.align_list[0].contig_pos2))
                fp.write("{} {} {} {}\n".format(
                    mis.align_list[1].ref_pos1, mis.align_list[1].ref_pos2,
                    mis.align_list[1].contig_pos1, mis.align_list[1].contig_pos2))

            print("Done")
            print()


def read_simpl_output(mis_class, args):
    """ Read output data from simplification stage """
    print("Reading output data from simplification stage...")
    try:
        fp = open(args.spades + SIMPLIFICATION_OUTPUT_FILENAME, "r")
    except:
        print("ERROR: some error on reading {}".format(args.spades + SIMPLIFICATION_OUTPUT_FILENAME))
        exit(1)
    else:
        with fp:
            fp.readline()  # contigs path
            fp.readline()  # ref path

            # broken bone
            data = fp.readline().split()
            cnt = int(data[0])
            for _ in range(cnt):
                bb = BrokenBoneSimplOutput()
                data = fp.readline().split()
                bb.is_ready = bool(int(data[0]))
                if bb.is_ready:
                    fp.readline()  # contig name
                    data = fp.readline().split()
                    bb.is_start = bool(int(data[0]))
                    if bb.is_start:
                        bb.is_start_broken = bool(int(data[1]))
                        if bb.is_start_broken:
                            bb.start_coverage = float(data[2])
                    data = fp.readline().split()
                    bb.is_end = bool(int(data[0]))
                    if bb.is_end:
                        bb.is_end_broken = bool(int(data[1]))
                        if bb.is_end_broken:
                            bb.end_coverage = float(data[2])
                mis_class.spades_broken_bone.append(bb)

            print("Done")
            print()


def print_results(mis_class, mis_list):
    """ Print results """
    print("================================================")
    print("RESULTS")
    print("{} misassemblies have been classified in this way:".format(len(mis_list)))
    print()

    cnt = len(mis_class.broken_bone)
    print("Broken bone:  {}\n".format(cnt))
    for i in range(cnt):
        print("{} of {}".format(i+1, cnt))
        print(mis_class.broken_bone[i])
        bb = mis_class.spades_broken_bone[i]
        if not bb.is_ready:
            print("Error has happened while processing the broken bone")
            continue
        if not bb.is_start:
            print("Start of the bone is correct")
        else:
            if bb.is_start_broken:
                print("Start of the bone has been broken. Ref2 incoming edge has been deleted due to low coverage {}".format(bb.start_coverage))
            else:
                print("Start of the bone has been broken. Unknown reason")
        if not bb.is_end:
            print("End of the bone is correct")
        else:
            if bb.is_end_broken:
                print("End of the bone has been broken. Ref1 outgoing edge has been deleted due to low coverage {}".format(bb.end_coverage))
            else:
                print("End of the bone has been broken. Unknown reason")
        print()
    print("------------------------------------------------\n")

    print("Ignored:  {}\n".format(len(mis_class.ignored)))
    print("------------------------------------------------\n")
    print("Unknown:  {}\n".format(len(mis_class.unknown)))



def parse_contigs(contigs_path):
    """ Parse contigs with extensive misassemblies """
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


def get_args():
    """ Parse and validate input arguments """
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

write_dist_est_input(mis_class, args)

print("Preparing SPAdes to start...")
print("Sorry, we are not able to set up SPAdes configs properly")
print("Please, set up 'configs/debruijn/config.info' manually")
print("Make sure 'entry_point distance_estimation' is selected")
print()
while input("Type 'ok' when configs are ready\n") != "ok":
    pass

print()
print("Starting SPAdes at the distance estimation stage")
print("====================================================")
subprocess.call([args.spades + "/run", "rd"])
print("====================================================")

print("Please, set up configs for 'entry_point simplification' to be selected")
print()
while input("Type 'ok' when configs are ready\n") != "ok":
    pass

print()
print("Starting SPAdes at the simplification stage")
print("====================================================")
subprocess.call([args.spades + "/run", "rd"])
print("====================================================")


read_simpl_output(mis_class, args)

print_results(mis_class, mis_list)

# delete temporary files
os.remove(DISTANCE_ESTIMATION_INPUT_FILENAME)
os.remove(DISTANCE_ESTIMATION_OUTPUT_FILENAME)
os.remove(SIMPLIFICATION_OUTPUT_FILENAME)