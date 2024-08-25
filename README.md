# snake-cma

snake-cma is Snakemake implementations for the Corrected Meta-Analysis (CMA) method.

It is developed and supported by team of scientists at the University of Edinburgh.

The method has the following features:
- It works for quantitative and binary traits, inheriting unique features of the extended method
- It is efficient; faster and using much less memory that existing GWAS software, suitable for very large sample sizes (> 1M)
- It can handle relatedness with large populations

There are currently two implementations:
- R+CMA, an implementation based on [REGENIE](https://rgcgithub.github.io/regenie)
- F+CMA: an implementation based on [fastGWA-GLMM](https://yanglab.westlake.edu.cn/software/gcta)


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
  -out OUT, --out Ousage: snake_cma.py  [<options> ...]

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
UT   output prefix
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
  --ncores NCORES       number of CPU cores to use
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
