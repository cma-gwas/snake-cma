# snake-cma

Corrected Meta-Analysis (CMA) is a resource-efficient method for genome-wide association studies (GWAS) on large populations.

snake-cma is CMA implementations on Snakemake.  

It is developed and supported by team of scientists at the University of Edinburgh.

It has the following features:
- It is efficient; faster and using much less memory that existing methods, suitable for very large sample sizes (> 1M)
- It works for quantitative and binary traits, inheriting unique features of the extended method
- It can handle relatedness with large populations

It currently supports and extends two existing methods:
- R+CMA, which uses [REGENIE](https://rgcgithub.github.io/regenie) for quantitative and binary traits
- F+CMA, which uses [fastGWA-GLMM](https://yanglab.westlake.edu.cn/software/gcta) for binary traits


## Quick start

- Download the source  
`$ git clone https://git.ecdf.ed.ac.uk/cma/snake-cma.git`
- Install dependencies  
`$ pip install snake-cma/requirements.txt`
- Run the main program  
`$ python snake-cma/src/snake_cma.py`

## Usage
```
usage: snake_cma.py  [<options> ...]

run CMA on Snakemake

options:
  -h, --help            show this help message and exit
  -bed BED, --bed BED   input bed prefix
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
  --ncores NCORES       number of cores to use
  -s SPLIT, --split SPLIT
                        number of sub-cohorts for splitting
  --bt                  binary trait analysis, default=false
  --cma-inflate CMA_INFLATE
                        CMA's inflation factor, default=0.2

```

## License
Apache 2.0 license

## Contributors  
- Mulya Agung (magung@ed.ac.uk)
- Ismail Ozkaraca (ismail.ozkaraca@ed.ac.uk)
