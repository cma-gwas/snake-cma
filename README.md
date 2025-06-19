# snake-cma

Corrected Meta-Analysis (CMA) is a resource-efficient method for genome-wide association studies (GWAS) on large populations.

**snake-cma** is the CMA implementation on Snakemake.  

It has the following features:
- It is efficient; faster and using much less memory than existing methods, suitable for very large sample sizes (> 1M)
- It works for quantitative and binary traits, inheriting unique features of the extended method
- It can handle relatedness with large populations

It currently supports and extends two existing methods:
- R+CMA, which uses [REGENIE](https://rgcgithub.github.io/regenie) for quantitative and binary traits
- F+CMA, which uses [fastGWA-GLMM](https://yanglab.westlake.edu.cn/software/gcta) for binary traits

## Quick start

- Download the source  
`$ git clone https://github.com/cma-gwas/snake-cma.git`
- Install dependencies with one of the following methods.  
  1. Create a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/installing-with-conda.html) environment using the provided environment.yml      
  `$ conda env create -f snake-cma/environment.yml`  
  `$ conda activate snake-cma`
  2. Use [Pip](https://packaging.python.org/en/latest/tutorials/installing-packages/) to install dependencies (**Python≥3.11** are supported)  
  `$ pip install -r snake-cma/requirements.txt`  
  > CMA also needs [PLINK](https://www.cog-genomics.org/plink/2.0/) and [REGENIE](https://rgcgithub.github.io/regenie/install/)/[GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Download) program binaries. See --plink, -rcma, and -fcma flags.
- Run the main program:  
  `$ cd snake-cma/src`  
  `$ python snake_cma.py`

## Usage
```
usage: snake_cma.py  [<options> ...]

Run CMA on Snakemake

options:
  -h, --help            show this help message and exit
  -bed BED, --bed BED   input genetic data, PLINK bed/bim/fam prefix
  -out OUT, --out OUT   output prefix
  -phen PHENOFILE, --phenoFile PHENOFILE
                        phenotype file
  -covar COVARFILE, --covarFile COVARFILE
                        covariate file
  -rcma REGENIE, --regenie REGENIE
                        REGENIE path for running R+CMA
  -fcma GCTA, --gcta GCTA
                        GCTA path for running F+CMA
  --sp-grm SP_GRM       sparse GRM path prefix for F+CMA
  --grm GRM             GRM path prefix for F+CMA, sparse GRM will be computed
  --plink PLINK         PLINK2 path for running F+CMA
  -s SPLIT, --split SPLIT
                        number of sub-cohorts for splitting
  --bt                  binary trait analysis, default=false
  --cma-inflate CMA_INFLATE
                        CMA inflation factor, default=0.2
  --ncores NCORES       number of cores to use
  --s-parallel S_PARALLEL
                        number of sub-cohorts in parallel, default=--split

```

## Citing
Mustafa İsmail Özkaraca, Mulya Agung, Pau Navarro, Albert Tenesa, Divide and conquer approach for genome-wide association studies, Genetics, Volume 229, Issue 4, April 2025, iyaf019, https://doi.org/10.1093/genetics/iyaf019

## License
CMA is developed at the University of Edinburgh and licensed under Apache License 2.0.

## Contributors  
- Mulya Agung (agung@hpc.is.tohoku.ac.jp)
- Ismail Ozkaraca (ismail.ozkaraca@ed.ac.uk)
