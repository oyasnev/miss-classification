import argparse
import os

__author__ = 'oyasnev'

QUAST_CONTIGS = "/contigs_reports/contigs_report_contigs.stdout"


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
