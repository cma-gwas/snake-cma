import os
import sys
from snakemake.logging import logger

# Mandatory environment vars
envvars:
    "PLINK_PATH",
    "REGENIE_PATH",
    "CMA_PATH",
    "SPLITTER_PATH",
    "WORKDIR",
    "INPUT_PREFIX",
    "PHENO_FILE",

## TOOLS ##
PLINK2 = os.getenv("PLINK_PATH")
REGENIE = os.getenv("REGENIE_PATH")
#META_PROG = PARENT_DIR + "/git-uoe/meta-analysis/meta.analysis.py"
META_PROG = os.getenv("CMA_PATH")
#SPLIT_PROG = PARENT_DIR + "/git-uoe/spark-meta/snakemake/cohort_splitter.py"
SPLIT_PROG = os.getenv("SPLITTER_PATH")
#COUNTER_SERVER_PROG = PARENT_DIR + "/git-uoe/spark-meta/snakemake/generic/counter_server.py"

# MANDATORY PARAMETERS
BED_PREFIX = os.getenv("INPUT_PREFIX")
OUT_DIR = os.getenv("WORKDIR")
PHENO = os.getenv("PHENO_FILE")

# OPTIONAL PARAMETERS
S1_ARGS = os.getenv("STEP1_ARGS")
S2_ARGS = os.getenv("STEP2_ARGS")
S1_SNPLIST = os.getenv("STEP1_SNPLIST")
S1_PRED_FILE_PREFIX = os.getenv("STEP1_PRED_FILE_PREFIX")
OVERLAP_CASE_SAMPLES = os.getenv("OVERLAP_CASE_SAMPLES", '')

S1_PRED_FILE_SUFFIX = "_pred.list"
extra_args = ""
split_extra_args = ""

COVAR_FILE = os.getenv("COVAR_FILE")
if COVAR_FILE:
    extra_args += " --covarFile {}".format(COVAR_FILE)

# Regenie only support one --keep parameter.
# So, ensure this will be accounted when creating the split id file.
# See META_ID_FILE
#ID_LIST = os.getenv("ID_LIST")
#if ID_LIST:
#    extra_args += " --keep {}".format(ID_LIST)

BINARY_TRAITS = os.getenv("BINARY_TRAITS", '')
if BINARY_TRAITS.lower() == "true":
    extra_args += " --bt"
    
    split_extra_args += " --binary-traits {}".format(PHENO)
    if OVERLAP_CASE_SAMPLES.lower() == "true":
        split_extra_args += " --overlap-case-samples"

s1_extra_args = ""
if S1_ARGS:
    s1_extra_args += S1_ARGS

s2_extra_args = ""
if S2_ARGS:
    s2_extra_args += S2_ARGS

SNP_LIST = os.getenv("SNP_LIST")
if SNP_LIST:
    s2_extra_args += " --extract {}".format(SNP_LIST)

if S1_SNPLIST:
    s1_extra_args += " --extract {}".format(S1_SNPLIST)

NUM_SPLITS = os.getenv("META_N_SPLITS")
if NUM_SPLITS:
  NUM_SPLITS = int(NUM_SPLITS)
else:
  NUM_SPLITS = 4

BISECT_RATIO = os.getenv("META_BISECT_RATIO")
if BISECT_RATIO and 0 < float(BISECT_RATIO) < 1.0:
    #print("NUM_SPLITS is ignored because META_BISECT_RATIO is set")
    split_extra_args += " --bisect-ratio {}".format(BISECT_RATIO)
    NUM_SPLITS = 2

SORTED_ID_LIST = os.getenv("SORTED_ID_LIST")
if SORTED_ID_LIST:
    split_extra_args += " --sorted-id-list {}".format(SORTED_ID_LIST)

PHEN_PREFIX = os.getenv("META_PHEN_LABEL_PREFIX")
if PHEN_PREFIX is None:
  PHEN_PREFIX = "Y"

PHEN_COLS = os.getenv("META_PHEN_COLS")
if PHEN_COLS: 
    PHEN_LIST = PHEN_COLS.split(',')
else:
    PHEN_LIST = list()

RAND_SAMPLES = os.getenv("META_RAND_SAMPLES")
if str(RAND_SAMPLES).lower() == "true":
    SPLIT_PROG += " -r"

RAND_SEED = os.getenv("META_RAND_SEED")
if RAND_SEED:
    RAND_SEED = int(RAND_SEED)
else:
    RAND_SEED = 123456789

QC_ARGS = os.getenv("META_QC_ARGS")
if QC_ARGS:
    QC_ARGS = QC_ARGS.strip()
else:
    QC_ARGS = ""

ID_LIST = os.getenv("ID_LIST")
if ID_LIST:
    QC_ID_FILE = ID_LIST
else:
    QC_ID_FILE = os.getenv("META_ID_FILE")

if QC_ID_FILE:
    QC_ARGS += " --keep {}".format(QC_ID_FILE)

## Input files ##
BED = BED_PREFIX + ".bed"

# Intermediary results
S1_OUT_PREFIX = OUT_DIR + "/step1_"
S2_OUT_PREFIX = OUT_DIR + "/step2_"
META_OUT_PREFIX = OUT_DIR + "/meta_"
QC_PREFIX = OUT_DIR + "/qc_pass"

# Phenotype names in the phenotype file
import pandas as pd

phen_headers = pd.read_csv(PHENO, index_col=False, nrows=0, sep="\s+").columns.tolist()
#NUM_PHENO = len(phen_headers) - 2
PHEN_IDS = phen_headers[2:]

if PHEN_LIST:
    PHEN_IDS = [p for p in PHEN_LIST if p in PHEN_IDS]

if not PHEN_IDS:
    print("Error: no phenotype names found in the phenotype file")
    sys.exit(-1)
else:
    extra_args += " --phenoColList {}".format(','.join(PHEN_IDS))

NUM_PHENO = len(PHEN_IDS)
#PHEN_IDS = [PHEN_PREFIX + str(x + 1) for x in range(NUM_PHENO)]

SUBCOHORT_IDS = [str(x+1) for x in range(NUM_SPLITS)]
# Split command pads zeroes
SPLIT_IDS = [s.zfill(2) for s in SUBCOHORT_IDS]
# Splits' filenames, results of the split command
SAMPLE_IDS = expand(QC_PREFIX + '.id.split.{subcohort_id}', subcohort_id=SUBCOHORT_IDS)

# Meta results
META_OUTPUTS = expand(META_OUT_PREFIX + '{phen_id}.meta', phen_id=PHEN_IDS)

# Count sample size for each subcohort
#import math
#import subprocess

#if ID_FILE:
#    with open(ID_FILE, 'r') as fp:
#        for count, line in enumerate(fp):
#            pass
#    NUM_SAMPLES = count + 1
#else:
#    NUM_SAMPLES =  int(subprocess.check_output(["wc", "-l", BED_PREFIX + ".fam"]).split()[0])

#NUM_LINES_EACH = math.ceil(NUM_SAMPLES / NUM_SPLITS)

# Schedule threads
NUM_PARALLEL = os.getenv("META_NUM_PARALLEL")
if NUM_PARALLEL:
    NUM_PARALLEL = int(NUM_PARALLEL.strip())
    if NUM_PARALLEL <= 0 or NUM_PARALLEL > NUM_SPLITS:
        NUM_PARALLEL = NUM_SPLITS
else:
    NUM_PARALLEL = NUM_SPLITS

step_nthreads = max(1, workflow.cores // NUM_PARALLEL)
meta_nthreads = max(1, workflow.cores // NUM_PHENO)

# Respect CPU affinities given by the system
cpuset = list(os.sched_getaffinity(0))
n_sub, rem = divmod(len(cpuset), NUM_PARALLEL)
sub_cpus = [cpuset[i:i+n_sub] for i in range(0, len(cpuset), n_sub)]
if rem > 0:
    sub_cpus[NUM_PARALLEL-1] += sub_cpus[-1]
    sub_cpus = sub_cpus[:NUM_PARALLEL]

# Connect to task counter server
import zmq
context = zmq.Context()
print("Connecting to counter server…")
socket = context.socket(zmq.REQ)
COUNTER_SERVER_PORT = os.getenv("COUNTER_SERVER_PORT", "5555")
socket.connect("tcp://localhost:{}".format(COUNTER_SERVER_PORT))

def counter_next(job_id):
    socket.send(bytes(job_id, encoding="utf-8"))
    b_value = socket.recv()
    return int.from_bytes(b_value, "big")

def counter_quit():
    socket.send(bytes("quit", encoding="utf-8"))

def get_sub_cpulist(task_id):
    job_counter = counter_next(task_id)
    job_cpuset = sub_cpus[job_counter % NUM_PARALLEL]
    logger.info("Task-{} cpuset: {}".format(task_id, job_cpuset))
    return ",".join([str(x) for x in job_cpuset])

onstart:
    #shell("python {COUNTER_SERVER_PROG} &")
    logger.info("Number of splits = {} (parallel={})".format(NUM_SPLITS, NUM_PARALLEL))

#onsuccess:
    #counter_next("quit")
#    counter_quit()

#onerror:
#    counter_quit()
    

rule all:
    input: META_OUTPUTS


rule qc:
    input: BED
    output: QC_PREFIX + ".id", QC_PREFIX + ".snplist"
    threads: workflow.cores
    shell: "{PLINK2} --bfile {BED_PREFIX} --write-snplist --write-samples --no-id-header --out {QC_PREFIX} --threads {threads} {QC_ARGS}"


rule split_sample_ids:
    input: QC_PREFIX + ".id"
    params:
        rand_seed=RAND_SEED
    output:
        SAMPLE_IDS
    threads: workflow.cores
    shell:
        "python {SPLIT_PROG} -s {QC_PREFIX}.id -o {QC_PREFIX}.id -n {NUM_SPLITS} --rand_seed {params.rand_seed} --threads {threads} {split_extra_args}"


rule step_one:
    input:
        bed=BED,
        pheno=PHENO,
        id=QC_PREFIX + ".id.split.{subcohort_id}"
    output:
        S1_OUT_PREFIX + '{subcohort_id}_pred.list'
    threads: step_nthreads
    run:
        # Reuse Step 1 predictor when provided
        pred_exists = False
        if S1_PRED_FILE_PREFIX:
            split_pred_file = S1_PRED_FILE_PREFIX + wildcards.subcohort_id + S1_PRED_FILE_SUFFIX
            if os.path.isfile(split_pred_file):
                print("Linking pred file {}->{}".format(split_pred_file, output[0]))
                shell("ln -s {split_pred_file} {output}")
                pred_exists = True

            #with open(split_pred_file, 'r') as f:
            #    loco_files = [line.rstrip().split(" ")[1] for line in f]
            #for loco in loco_files:
            #    shell("ln -s {loco} {output}")
            #    loco_basename, loco_dir = os.path.basename(loco), os.path.dirname(loco)
            #    os.makedirs(loco_dir,exist_ok=True)
        if not pred_exists:
            cpu_list = get_sub_cpulist(wildcards.subcohort_id)

            shell("""
                /usr/bin/time -v -o {S1_OUT_PREFIX}{wildcards.subcohort_id}.time taskset -c {cpu_list} {REGENIE} --step 1 --bed {BED_PREFIX} --bsize 1000 --phenoFile {input.pheno} --out {S1_OUT_PREFIX}{wildcards.subcohort_id} --keep {input.id} --threads {threads} {extra_args} {s1_extra_args}
            """)
    #shell:
    #    "/usr/bin/time -v -o {S1_OUT_PREFIX}{wildcards.subcohort_id}.time {REGENIE} --step 1 --bed {BED_PREFIX} --bsize 1000 --phenoFile {input.pheno} --out {S1_OUT_PREFIX}{wildcards.subcohort_id} --keep {QC_PREFIX}.id.split.{wildcards.subcohort_id} --threads {threads} {extra_args} {s1_extra_args}"
        #"/usr/bin/time -v -o {S1_OUT_PREFIX}{wildcards.subcohort_id}.time {REGENIE} --step 1 --bed {BED_PREFIX} --bsize 1000 --phenoFile {input.pheno} --out {S1_OUT_PREFIX}{wildcards.subcohort_id} --keep {QC_PREFIX}.id.split.{wildcards.subcohort_id} --extract {input.snp} --threads {threads} {extra_args} {s1_extra_args}"


rule step_two:
    input:
        bed=BED,
        pheno=PHENO,
        id=QC_PREFIX + '.id.split.{subcohort_id}',
        pred=S1_OUT_PREFIX + '{subcohort_id}_pred.list'
    output:
        #S2_OUT_PREFIX + '{subcohort_id}_{phen_id}.regenie'
        expand(S2_OUT_PREFIX + '{{subcohort_id}}_{phen_id}.regenie', phen_id=PHEN_IDS)
    threads: step_nthreads
    run:
        cpu_list = get_sub_cpulist(wildcards.subcohort_id)

        shell("""
            /usr/bin/time -v -o {S2_OUT_PREFIX}{wildcards.subcohort_id}.time taskset -c {cpu_list} {REGENIE} --step 2 --bed {BED_PREFIX} --bsize 2000 --phenoFile {input.pheno} --pred {input.pred} --out {S2_OUT_PREFIX}{wildcards.subcohort_id} --keep {input.id} --threads {threads} {extra_args} {s2_extra_args}
            """)
    #shell:
    #    "/usr/bin/time -v -o {S2_OUT_PREFIX}{wildcards.subcohort_id}.time {REGENIE} --step 2 --bed {BED_PREFIX} --bsize 2000 --phenoFile {input.pheno} --pred {input.pred} --out {S2_OUT_PREFIX}{wildcards.subcohort_id} --keep {QC_PREFIX}.id.split.{wildcards.subcohort_id} --extract {input.snp} --threads {threads} {extra_args} {s2_extra_args}"


rule meta_analysis:
    input:
        expand(S2_OUT_PREFIX + '{subcohort_id}_{{phen_id}}.regenie', subcohort_id=SUBCOHORT_IDS)
    output:
        META_OUT_PREFIX + '{phen_id}.meta'
    threads: meta_nthreads
    shell:
        "/usr/bin/time -v -o {META_OUT_PREFIX}{wildcards.phen_id}.time python {META_PROG} --meta_summary {S2_OUT_PREFIX}%s_{wildcards.phen_id}.regenie --split {NUM_SPLITS} --inflate 0.2 --out {META_OUT_PREFIX}{wildcards.phen_id}"

