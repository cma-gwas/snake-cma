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




