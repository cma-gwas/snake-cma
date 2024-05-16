import argparse
import time
from functions import *
from import_config import load_config

# Argparse
parser = argparse.ArgumentParser(prog="META-ANALYSIS", description="Meta-Analysis Software", epilog= "For more information please see https://git.ecdf.ed.ac.uk/iozkarac/meta-analysis.git")
parser.add_argument("--meta_summary", help="Example Usage: --meta_summary /home/Ozkaraca/summary%%s.csv", metavar="<files_path_prefix>")
parser.add_argument("--meta_no_correction", help="Example Usage: --meta_no_correction /home/Ozkaraca/summary%%s.csv", metavar="<files_path_prefix>")
parser.add_argument("--split", help="Example Usage: --split 16", metavar="<split_index>")
parser.add_argument("--out", help="output help", default="summary.stats", metavar="<output file name>")
parser.add_argument("--inflate", help="Example Usage: --inflate 0.5", metavar="<prior heritability estimate of the trait>")
parser.add_argument("--binary_traits", action="store_true", help="Binary traits analysis, default is quantitative")
#parser.add_argument("--is_binary", help="Example Usage: --is_binary 1 (trait is a disease), or is_binary 0 (trait is continuous)", metavar="<0 or 1>", default="0")
parser.add_argument("--fastgwa", action="store_true", help="Use fastGWA's summary file format, default is REGENIE's")

initial_message='''
A Meta-Analysis Tool..
A manual is available here:
https://git.ecdf.ed.ac.uk/iozkarac/meta-analysis.git
__________________________________
'''
initial_log='''
META-ANALYSIS v.0.2

Analysis Start Time: %s

Working Directory: %s

Input Variables:

'''


if __name__ == "__main__":
    start_time=time.time()
    print(initial_message)
    args = parser.parse_args()
    initial_logging(initial_log, vars(parser.parse_args()), vars(parser.parse_args([])), args.out + ".log", args.split, args.meta_summary, args.meta_no_correction, args.out)

    if args.fastgwa:
        import_conf = load_config("fastgwa")
    else:
        import_conf = load_config("regenie")

#    if args.meta_summary:
#        print('\nPerforming fixed-effect meta analysis with correction using summary statistics...\nSplit index is: %s' % args.split)
#        X = meta_summary(args.meta_summary, int(args.split))
#        X.to_csv(args.out + ".meta", sep=" ")

    if args.meta_summary:
        print('\nPerforming fixed-effect meta analysis with correction using summary statistics...\nSplit index is: %s' % args.split)
        if args.inflate:
            X = meta_summary_inflated(args.meta_summary, int(args.split), float(args.inflate), args.binary_traits,
                                      import_conf)
            X.to_csv(args.out + ".meta", sep=" ")
        else:
            X = meta_summary(args.meta_summary, int(args.split), import_conf)
            X.to_csv(args.out + ".meta", sep=" ")

    if args.meta_no_correction :
        print('\nPerforming fixed-effect meta analysis without any correction...\nSplit index is: %s' % args.split)
        X = meta_nocorrection(args.meta_no_correction, int(args.split), import_conf)
        X.to_csv(args.out + ".meta", sep=" ")

    with open(args.out + ".log", "a") as f:
        f.write("\nAnalysis End Time: %s" % datetime.datetime.now().strftime("%A - %d %B %Y - %H:%M:%S"))
        f.write("\nProgram halted in %s seconds...\n" % (round(time.time() - start_time,2)))


