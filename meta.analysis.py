import argparse
import os
import time
import datetime
from functions import *

# Argparse
parser = argparse.ArgumentParser(prog="META-ANALYSIS", description="Meta-Analysis Software", epilog= "For more information please see https://git.ecdf.ed.ac.uk/iozkarac/meta-analysis.git")
parser.add_argument("--meta_summary", help="Example Usage: --meta_summary /home/Ozkaraca/summary%%s.csv", metavar="<files_path_prefix>")
parser.add_argument("--meta_no_correction", help="Example Usage: --meta_no_correction /home/Ozkaraca/summary%%s.csv", metavar="<files_path_prefix>")
parser.add_argument("--split", help="Example Usage: --split 16", metavar="<split_index>")
parser.add_argument("--out", help="output help", default="summary.stats", metavar="<output file name>")

initial_message='''
A Meta-Analysis Tool..
A manual is available here:
https://git.ecdf.ed.ac.uk/iozkarac/meta-analysis.git
__________________________________
'''
initial_log='''
META-ANALYSIS v.0.1

Analysis Start Time: %s

Working Directory: %s

Input Variables:

'''


if __name__ == "__main__":
    start_time=time.time()
    print(initial_message)
    args = parser.parse_args()
    initial_logging(initial_log, vars(parser.parse_args()), vars(parser.parse_args([])), args.out + ".log")

    if args.meta_summary:
        print('\nPerforming fixed-effect meta analysis with correction using summary statistics...\nSplit index is: %s' % args.split)
        X = meta_summary(args.meta_summary, int(args.split))
        X.to_csv(args.out + ".meta", sep=" ")

    if args.meta_no_correction :
        print('\nPerforming fixed-effect meta analysis without any correction...\nSplit index is: %s' % args.split)
        X = meta_nocorrection(args.meta_no_correction, int(args.split))
        X.to_csv(args.out + ".meta", sep=" ")

    with open(args.out + ".log", "a") as f:
        f.write("\nAnalysis End Time: %s" % datetime.datetime.now().strftime("%A - %d %B %Y - %H:%M:%S"))
        f.write("\nProgram halted in %s seconds...\n" % (round(time.time() - start_time,2)))


