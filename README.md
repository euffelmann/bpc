# BPC

The Bayesian polygenic score Probability Conversion (BPC) approach transforms Bayesian polygenic scores (PGSs) into predicted disorder probabilities. The BPC approach is described in the following preprint: https://www.medrxiv.org/content/10.1101/2024.01.12.24301157v1.

### Getting Started

The following files and tools will need to be downloaded:

-   The `bpc.R` and `h2o_to_h2l.R` files. Alternatively you can clone this whole repository using the following git command: `git clone https://github.com/euffelmann/bpc.git`. `toy_data.profile` contains simulated Polygenic Scores (PGS) in plink format. 

-   1000 Genomes population reference files (from e.g. <https://ctg.cncr.nl/software/magma>)

-   PRScs (<https://github.com/getian107/PRScs>) or GCTB for SBayesR (<https://cnsgenomics.com/software/gctb/#Download>)

-   plink1.9 to compute PGSs (<https://www.cog-genomics.org/plink/>)

### Applying the BPC function

The BPC function takes four inputs:

-   **pgs_liab**: A single individual's PGS (or a vector of PGSs) on the liability scale based on the posterior mean betas from a Bayesian PGS method. While in theory any Bayesian PGS method that is well-calibrated for continuous traits can be used, we have only evaluated PRScs and SBayesR.

    -   PRScs: Use the effective sample size (`neff = 4 / ((1 / n_cases) + (1 / n_controls))`) of the GWAS training sample as input to compute posterior mean betas. If the GWAS was based on multiple cohorts, use the sum of all cohorts' effective sample sizes (see Grotzinger et al. (2023) Biological Psychiatry, https://doi.org/10.1016/j.biopsych.2022.05.029). After calculating the PGSs with plink1.9, transform them from the observed scale with 50% case ascertainment (because neff is used as input to PRScs) to the liability scale. Example code:

        ```         
        ## calculate PGSs with plink based on default PRScs output and 
        ## 1000 Genomes allele frequencies to shift all scores to mean zero
        ## command line (example code, the output is provided in toy_data/):
        plink1.9 \ 
          --bfile <filename> \ 
          --read-freq g1000_eur.frq \ # can be computed with 'plink1.9 --freq'
          --score <PRScs output filename> 2 4 6 sum center \
          --out toy_data/toy_testdata

        ## transform PGSs from the observed to the liability scale
        ## R code:
        library(data.table)
        source("h2o_to_h2l.R")
        pgs_obs <- read.table("toy_data/toy_testdata.profile", header = TRUE)
        # for a disorder with a population prevalence of 0.01. adjust K accordingly
        pgs_liab <- pgs_obs$SCORESUM * sqrt(h2o_to_h2l(K = 0.01))
        ```

    -   SBayesR: To achieve well-calibrated PGSs, the GWAS betas and standard errors need to be transformed to the standardized observed scale with 50% case ascertainment first before SBayesR can be used.

        ```         
        ## R code:         
        sumstats$b5050 <- sumstats$b / (sumstats$se * sqrt(sumstats$neff))
        sumstats$se <- 1 / sqrt(sumstats$neff)
        ```

        After the GWAS betas and standard errors have been transformed in this way, you can use SBayesR with neff as input and rescale the PGS to the liability scale in the same way as for PRScs. Example code:

        ```         
        ## compute PGSs with plink based on default SBayesR output and 
        ## 1000 Genomes allele frequencies to shift all scores to mean zero
        ## command line (example code, the output is provided in toy_data/):
        plink1.9 \ 
          --bfile <filename> \ 
          --read-freq g1000_eur.frq \ # can be computed with 'plink1.9 --freq'
          --score <SBayesR output filename> 2 5 8 sum center \
          --out toy_data/toy_testdata

        ## transform PGSs from the observed to the liability scale
        ## R code:
        library(data.table)
        source("h2o_to_h2l.R")
        pgs_obs <- read.table("toy_data/toy_testdata.profile", header = TRUE)
        pgs_liab <- pgs_obs$SCORESUM * sqrt(h2o_to_h2l(K = 0.01))
        ```

-   **K**: The population prevalence of the disorder of interest (e.g. 0.01 for schizophrenia)

-   **P**: The prior disorder probability. This can be interpreted as the proportion of cases in a (hypothetical) testing sample.

-   **r2l**: The explained variance (*R*<sup>2</sup>) on the liability scale. In the BPC approach, the variance of a well-calibrated PGS in a population reference sample (e.g. 1000 Genomes) is used to estimate *R*<sup>2</sup><sub>liability</sub>. Simply use **the same** posterior mean betas as used above to calculate PGSs in 1000 Genomes, then transform the PGSs from the observed to the liability scale, and compute the variance of the PGSs to estimate *R*<sup>2</sup><sub>liability</sub>. Example code with K = 0.01:

    ```         
    ## transform PGSs from the observed to the liability scale in 1000 Genomes
    ## R code:
    pgs_ref_obs <- read.table("toy_data/toy_refdata.profile", header = TRUE)
    pgs_ref_liab <- pgs_ref_obs$SCORESUM * sqrt(h2o_to_h2l(K = 0.01))

    r2l = var(pgs_ref_liab)
    ```

After preparing the input, the BPC function can be applied as follows to compute the predicted disorder probabilities (for a population prevalence of 0.01 and a prior disorder probability of 0.5):

```         
## R code:
source("bpc.R")
pred_prob <- bpc(pgs_liab = pgs_liab, K = 0.01, prior = 0.5, r2l = r2l)
```

The output is a vector of predicted disorder probabilities for the individual(s) for whom the PGSs were computed. See `toy_data/toy_pred_prob.txt` for the expected output (PHENO = case control status, pred_prob = predicted disorder probabilities).


We note that the prior for a random individual from the full population would be equal to the population prevalence, e.g. 0.01 for schizophrenia, but the prior for a help-seeking individual in a clinical setting will be substantially higher. Furthermore, the BPC approach has only been tested in individuals of European ancestry. Ancestry mismatches between the GWAS training sample and the individual for whom the BPC is computed may negatively impact its calibration.

#### Depencies

The R and package versions used can be found in the `renv.lock` file.