import argparse
import os
import time
import datetime

# Argparse
parser = argparse.ArgumentParser(prog="META-ANALYSIS", description="Meta-Analysis Software", epilog= "For more information please see https://git.ecdf.ed.ac.uk/iozkarac/meta-analysis.git")
parser.add_argument("--meta_summary", help="Example Usage: --meta_summary /home/Ozkaraca/summary%%s.csv", metavar="<files_path_prefix>")
parser.add_argument("--meta_no_correction", help="Example Usage: --meta_no_correction /home/Ozkaraca/summary%%s.csv", metavar="<files_path_prefix>")
parser.add_argument("--split", help="Example Usage: --split 16", metavar="<split_index>")
parser.add_argument("--out", help="output help", default="ace_file", metavar="<output file name>")

initial_message='''
A Meta-Analysis Tool..

A manual is available here:
https://git.ecdf.ed.ac.uk/iozkarac/meta-analysis.git
__________________________________
'''

if __name__ == "__main__":
        from functions import *
        print(initial_message)
        args = parser.parse_args()
        input_vars = vars(parser.parse_args())
        with open(args.out + ".log", "w") as f:
                f.write("META-ANALYSIS v.0.0")
                f.write("\nAnalysis Start Time: %s" % datetime.datetime.now().strftime("%A - %d %B %Y - %H:%M:%S"))
                start_time=time.time()
                f.write("\nWorking Directory: %s" % str(os.getcwd()))
                f.write("\nInput Variables:\n")
                for key in input_vars:
                        if input_vars[key]!= None:
                                f.write("\t--%s %s\n" % (str(key), str(input_vars[key])))

        if args.meta_summary:
                print('\nPerforming fixed-effect meta analysis with correction using summary statistics...')
                split_index = int(args.split)
                print('\nSplit index is: %s' % args.split)
                path_prefix = args.meta_summary
                X = meta_summary(path_prefix, split_index)
                X.to_csv(args.out + ".out", sep=" ")
                del X

        if args.meta_no_correction :
                print('\nPerforming fixed-effect meta analysis without any correction...')
                split_index = int(args.split)
                print('\nSplit index is: %s' % args.split)
                path_prefix = args.meta_no_correction
                X = meta_nocorrection(path_prefix, split_index)
                X.to_csv(args.out + ".out", sep=" ")
                del X

        with open(args.out + ".log", "a") as f:
                f.write("\nAnalysis End Time: %s" % datetime.datetime.now().strftime("%A - %d %B %Y - %H:%M:%S"))
                f.write("\nProgram halted in %s seconds...\n" % (round(time.time() - start_time,2)))
        del input_vars


