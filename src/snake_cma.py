import sys
import os
import argparse
import subprocess

import psutil
from snakemake import api
from snakemake import settings
from pathlib import Path

this_dir = os.path.abspath(os.path.dirname(__file__))
rcma_snakefile = os.path.join(this_dir, "snakefiles", "rcma.Snakefile")
fcma_snakefile = os.path.join(this_dir, "snakefiles", "fcma.Snakefile")
cma_prog = os.path.join(this_dir, "meta_analysis", "meta.analysis.py")
splitter_prog = os.path.join(this_dir, "cohort_splitter.py")
counter_server = os.path.join(this_dir, "sched", "counter_server.py")

def main():
    parser = argparse.ArgumentParser(description='run CMA on Snakemake', usage='''snake_cma.py  [<options> ...]''')

    parser.add_argument('-bed', '--bed', type=str, help='input bed prefix')
    parser.add_argument('-out', '--out', type=str, help='output prefix', required=True)
    parser.add_argument('-phen', '--phenoFile', type=str, help='phenotype file', required=True)
    parser.add_argument('-covar', '--covarFile', type=str, help='covariate file')
    parser.add_argument('-rcma', '--regenie', type=str, help='REGENIE path for running R+CMA')
    parser.add_argument('-fcma', '--gcta', type=str, help='GCTA path for running F+CMA')
    parser.add_argument('--sp-grm', type=str, help='sparse GRM path prefix for F+CMA')
    parser.add_argument('--grm', type=str, help='GRM path prefix for F+CMA, sparse GRM will be computed')
    parser.add_argument('--plink', type=str, help='PLINK2 path for running F+CMA', required=True)
    parser.add_argument('-s', '--split', type=int, help='number of sub-cohorts for splitting', default=4)
    parser.add_argument( '--bt', action="store_true", help='binary trait analysis, default=false')
    parser.add_argument('--cma-inflate', type=float, help="CMA inflation factor, default=0.2")
    parser.add_argument('--ncores', type=int, help='number of cores to use',
                        default=psutil.cpu_count(logical=False))
    parser.add_argument('--s-parallel', type=int,
                        help='number of sub-cohorts in parallel, default=--split')

    args = parser.parse_args()

    # Program path parameters
    #os.environ["COUNTER_SERVER_PATH"] = counter_server
    os.environ["SPLITTER_PATH"] = splitter_prog
    os.environ["CMA_PATH"] = cma_prog
    os.environ["PLINK_PATH"] = args.plink
    os.environ["WORKDIR"] = args.out

    # Analysis parameters
    os.environ["INPUT_PREFIX"] = args.bed
    os.environ["PHENO_FILE"] = args.phenoFile
    os.environ["META_N_SPLITS"] = str(args.split)

    if args.bt:
        os.environ["BINARY_TRAITS"] = "true"

    if args.regenie:
        os.environ["REGENIE_PATH"] = args.regenie
        snakefile = rcma_snakefile
    elif args.gcta and (args.sp_grm or args.grm) and args.bt:
        os.environ["GCTA_PATH"] = args.gcta
        if args.sp_grm:
            os.environ["SP_GRM_FILE_PREFIX"] = args.sp_grm
        else:
            os.environ["GRM_FILE_PREFIX"] = args.grm
        snakefile = fcma_snakefile
    else:
        sys.stderr.write("Error: either -rcma or -fcma must be specified, and -fcma needs --grm and --bt.\n")
        sys.exit(-1)

    if not os.path.exists(snakefile):
        sys.stderr.write("Error: cannot find {}\n".format(snakefile))
        sys.exit(-1)

    if args.covarFile:
        os.environ["COVAR_FILE"] = args.covarFile

    if args.cma_inflate:
        if args.cma_inflate > 0:
            os.environ["CMA_INFLATE"] = str(args.cma_inflate)
        else:
            sys.stderr.write("Error: --cma-inflate must be > 0.\n")
            sys.exit(-1)

    # Scheduling parameters
    if args.s_parallel:
        os.environ["META_NUM_PARALLEL"] = str(args.s_parallel)

    # Run workflow
    with api.SnakemakeApi(
        settings.OutputSettings(
            verbose=False,
            show_failed_logs=True,
        )
    ) as snakemake_api:
        rs = settings.ResourceSettings(
            cores=args.ncores,
        )
        workflow_api = snakemake_api.workflow(
            storage_settings=settings.StorageSettings(),
            resource_settings=rs,
            snakefile=Path(snakefile),
            workdir=Path(args.out),
        )
        dag_api = workflow_api.dag()

        # Run counter_server as a background process
        counter_server_pid = subprocess.Popen(["python", counter_server]).pid

        # Go on by calling methods of the dag api.
        try:
            dag_api.execute_workflow()
        except Exception as e:
            sys.stderr.write(e.__str__() + "\n")
            sys.exit(-1)
        finally:
            os.system("kill -9 {}".format(counter_server_pid))


if __name__ == "__main__":
    sys.exit(main())