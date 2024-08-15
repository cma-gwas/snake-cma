import os
from snakemake.logging import logger

# Mandatory environment vars
envvars:
    "PROG_PARENT_DIR",
    "META_WORKDIR",
    "INPUT_PREFIX",
    "PHENO_FILE",

## PARENT DIRECTORY ##
PARENT_DIR = os.getenv("PROG_PARENT_DIR")

## TOOLS ##
PLINK2 = PARENT_DIR + "/prog/plink2"
GCTA = PARENT_DIR + "/prog/gcta64"
META_PROG = PARENT_DIR + "/git-uoe/meta-analysis/meta.analysis.py"
SPLIT_PROG = PARENT_DIR + "/git-uoe/spark-meta/snakemake/cohort_splitter.py"
COUNTER_SERVER_PROG = PARENT_DIR + "/git-uoe/spark-meta/snakemake/generic/counter_server.py"

# MANDATORY PARAMETERS
BED_PREFIX = os.getenv("INPUT_PREFIX")
OUT_DIR = os.getenv("META_WORKDIR")
PHENO = os.getenv("PHENO_FILE")

# OPTIONAL PARAMETERS
COVAR_FILE = os.getenv("COVAR_FILE")
OVERLAP_CASE_SAMPLES = os.getenv("OVERLAP_CASE_SAMPLES", '')
S1_ARGS = os.getenv("STEP1_ARGS")
S2_ARGS = os.getenv("STEP2_ARGS")
GRM_FILE_PREFIX = os.getenv("GRM_FILE_PREFIX")
SP_GRM_FILE_PREFIX = os.getenv("SP_GRM_FILE_PREFIX")

SP_GRM_FILE_SUFFIX = ".grm.sp"
SP_GRM_ID_FILE_SUFFIX = ".grm.id"

split_extra_args = ""

s1_extra_args = ""
if S1_ARGS:
    s1_extra_args += S1_ARGS

s2_extra_args = ""
if S2_ARGS:
    s2_extra_args += S2_ARGS

QC_ARGS = os.getenv("META_QC_ARGS")
if QC_ARGS:
    QC_ARGS = QC_ARGS.strip()
else:
    QC_ARGS = ""

ID_FILE = os.getenv("ID_LIST")
if ID_FILE:
    QC_ARGS += " --keep {}".format(ID_FILE)
    #extra_args += " --keep {}".format(ID_LIST)

BINARY_TRAITS = os.getenv("BINARY_TRAITS", '')
if BINARY_TRAITS.lower() == "true":
    s2_extra_args += " --fastGWA-mlm-binary"
    
    split_extra_args += " --binary-traits {}".format(PHENO)
    if OVERLAP_CASE_SAMPLES.lower() == "true":
        split_extra_args += " --overlap-case-samples"

if COVAR_FILE:
    s2_extra_args += " --qcovar {}".format(COVAR_FILE)

SNP_LIST = os.getenv("SNP_LIST")
if SNP_LIST:
    s2_extra_args += " --extract {}".format(SNP_LIST)

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

RAND_SAMPLES = os.getenv("META_RAND_SAMPLES")
if str(RAND_SAMPLES).lower() == "true":
    SPLIT_PROG += " -r"

RAND_SEED = os.getenv("META_RAND_SEED")
if RAND_SEED:
    RAND_SEED = int(RAND_SEED)
else:
    RAND_SEED = 123456789

## Input files ##
BED = BED_PREFIX + ".bed"

# Intermediary results
S1_OUT_PREFIX = OUT_DIR + "/step1_"
S2_OUT_PREFIX = OUT_DIR + "/step2_"
META_OUT_PREFIX = OUT_DIR + "/meta_"
QC_PREFIX = OUT_DIR + "/qc_pass"

# Phenotype names in the phenotype file
import pandas as pd

phen_headers = pd.read_csv(PHENO,index_col=False, nrows=0, sep="\s+").columns.tolist()
NUM_PHENO = len(phen_headers) - 2

PHEN_IDS = [PHEN_PREFIX + str(x + 1) for x in range(NUM_PHENO)]

SUBCOHORT_IDS = [str(x+1) for x in range(NUM_SPLITS)]
# Split command pads zeroes
SPLIT_IDS = [s.zfill(2) for s in SUBCOHORT_IDS]
# Splits' filenames, results of the split command
SAMPLE_IDS = expand(QC_PREFIX + '.id.split.{subcohort_id}', subcohort_id=SUBCOHORT_IDS)

# Meta results
META_OUTPUTS = expand(META_OUT_PREFIX + '{phen_id}.meta', phen_id=PHEN_IDS)

# Count sample size for each subcohort
import math
import subprocess

if ID_FILE:
    with open(ID_FILE, 'r') as fp:
        for count, line in enumerate(fp):
            pass
    NUM_SAMPLES = count + 1
else:
    NUM_SAMPLES =  int(subprocess.check_output(["wc", "-l", BED_PREFIX + ".fam"]).split()[0])

NUM_LINES_EACH = math.ceil(NUM_SAMPLES / NUM_SPLITS)

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
socket.connect("tcp://localhost:5555")

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
        id=QC_PREFIX + ".id.split.{subcohort_id}",
        grm=GRM_FILE_PREFIX + ".grm.bin"
    output:
        sp_grm=S1_OUT_PREFIX + '{subcohort_id}' + SP_GRM_FILE_SUFFIX,
        sp_grm_id=S1_OUT_PREFIX + '{subcohort_id}' + SP_GRM_ID_FILE_SUFFIX
    params:
        sp_cutoff=0.05
    threads: step_nthreads
    run:
        # Reuse sp_grm when provided
        sp_grm_exists = False
        if SP_GRM_FILE_PREFIX:
            sp_grm_file = SP_GRM_FILE_PREFIX + wildcards.subcohort_id + SP_GRM_FILE_SUFFIX
            sp_grm_id_file = SP_GRM_FILE_PREFIX + wildcards.subcohort_id + SP_GRM_ID_FILE_SUFFIX
            if os.path.isfile(sp_grm_file) and os.path.isfile(sp_grm_id_file):
                print("Linking sp_grm file {}->{}".format(sp_grm_file, output[0]))
                shell("ln -s {sp_grm_file} {output.sp_grm}")
                shell("ln -s {sp_grm_id_file} {output.sp_grm_id}")
                sp_grm_exists = True

        if not sp_grm_exists:
            cpu_list = get_sub_cpulist(wildcards.subcohort_id)

            shell("""
                /usr/bin/time -v -o {S1_OUT_PREFIX}{wildcards.subcohort_id}.time taskset -c {cpu_list} {GCTA} --bfile {BED_PREFIX} --keep {input.id} --grm {GRM_FILE_PREFIX} --make-bK-sparse {params.sp_cutoff} --out {S1_OUT_PREFIX}{wildcards.subcohort_id} --threads {threads} {s1_extra_args}
            """)


rule step_two:
    input:
        bed=BED,
        pheno=PHENO,
        id=QC_PREFIX + '.id.split.{subcohort_id}',
        sp_grm=S1_OUT_PREFIX + '{subcohort_id}' + SP_GRM_FILE_SUFFIX
    output:
        expand(S2_OUT_PREFIX + '{{subcohort_id}}_{phen_id}.fastGWA',phen_id=PHEN_IDS)
    params:
        sp_grm_prefix=S1_OUT_PREFIX + '{subcohort_id}'
    threads: step_nthreads
    run:
        cpu_list = get_sub_cpulist(wildcards.subcohort_id)

        for i in range(1, NUM_PHENO+1):
            phen_id = PHEN_IDS[i-1]
            shell("""
                /usr/bin/time -v -o {S2_OUT_PREFIX}{wildcards.subcohort_id}.time taskset -c {cpu_list} {GCTA} --bfile {BED_PREFIX} --pheno {input.pheno} --mpheno {i} --grm-sparse {params.sp_grm_prefix} --nofilter --out {S2_OUT_PREFIX}{wildcards.subcohort_id}_{phen_id} --keep {input.id} --threads {threads} {s2_extra_args}
                """)


rule meta_analysis:
    input:
        #expand(S2_OUT_PREFIX + '{subcohort_id}_{phen_id}.fastGWA',subcohort_id=SUBCOHORT_IDS,phen_id=PHEN_IDS)
        expand(S2_OUT_PREFIX + '{subcohort_id}_{{phen_id}}.fastGWA',subcohort_id=SUBCOHORT_IDS)
    output:
        META_OUT_PREFIX + '{phen_id}.meta'
    threads: meta_nthreads
    shell:
        "/usr/bin/time -v -o {META_OUT_PREFIX}{wildcards.phen_id}.time python {META_PROG} --meta_summary {S2_OUT_PREFIX}%s_{wildcards.phen_id}.fastGWA --split {NUM_SPLITS} --fastgwa --binary_traits --inflate 0.2 --out {META_OUT_PREFIX}{wildcards.phen_id}"

