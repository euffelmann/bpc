# BPC

The Bayesian polygenic score Probability Conversion (BPC) approach transforms Bayesian polygenic scores (PGSs) into predicted disorder probabilities. The BPC approach is described in the following paper: [TBA].

### Getting Started

The following files and tools will need to be downloaded:

-   The `bpc.R` and `h2o_to_h2l.R` files. Alternatively you can clone this whole repository using the following git command: `git clone https://github.com/euffelmann/bpc.git`

-   1000 Genomes population reference files (from e.g. <https://ctg.cncr.nl/software/magma>)

-   PRScs (<https://github.com/getian107/PRScs>) or GCTB for SBayesR (<https://cnsgenomics.com/software/gctb/#Download>)

-   plink1.9 to compute PGSs (<https://www.cog-genomics.org/plink/>)

### Applying the BPC function

The BPC function takes four inputs:

-   **pgs_liab**: A single individual's PGS (or a vector of PGSs) on the liability scale based on the posterior mean betas from a Bayesian PGS method. While in theory any Bayesian PGS method that is well-calibrated for continuous traits can be used, we have only evaluated PRScs and SBayesR.

    -   PRScs: Use the effective sample size (`neff = 4 / ((1 / n_cases) + (1 / n_controls`))) of the GWAS training sample as input to compute posterior mean betas. After calculating the PGSs with plink1.9, transform them from the observed scale with 50% case ascertainment to the liability scale. Example code:

        ```         
        ## calculate PGSs with plink based on default PRScs output and 
        ## 1000 Genomes allele frequencies to shift all scores to mean zero
        ## command line:
        plink1.9 \ 
          --bfile <filename> \ 
          --read-freq g1000_eur.frq \ # can be computed with 'plink1.9 --freq'
          --score <PRScs output filename> 2 4 6 sum center \
          --out pgs_obs

        ## transform PGSs from the observed to the liability scale
        ## R code:
        library(data.table)
        source(h2o_to_h2l.R)
        pgs_obs <- fread("pgs_obs.profile", header = TRUE)
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
        ## command line:
        plink1.9 \ 
          --bfile <filename> \ 
          --read-freq g1000_eur.frq \ # can be computed with 'plink1.9 --freq'
          --score <SBayesR output filename> 2 5 8 sum center \
          --out pgs_obs

        ## transform PGSs from observed scale to liability scale
        ## R code:
        library(data.table)
        source(h2o_to_h2l.R)
        pgs_obs <- fread("pgs_obs.profile", header = TRUE)
        pgs_liab <- pgs_obs$SCORESUM * sqrt(h2o_to_h2l(K = 0.01))
        ```

-   **K**: The population prevalence of the disorder of interest (e.g. 0.01 for schizophrenia)

-   **P**: The prior disorder probability. This can be interpreted as the proportion of cases in a (hypothetical) testing sample.

-   **r2l**: The explained variance (*R*<sup>2</sup>) on the liability scale. In the BPC approach, the variance of a well-calibrated PGS in a population reference sample (e.g. 1000 Genomes) is used to estimate *R*<sup>2</sup><sub>liability</sub>. Simply use **the same** posterior mean betas from PRScs or SBayesR to calculate PGSs in 1000 Genomes, then transform the PGSs from the observed to the liability scale, and compute the variance of the PGSs to estimate *R*<sup>2</sup><sub>liability</sub>. Example code with K = 0.01:

    ```         
    ## transform PGSs from the observed to the liability scale in 1000 Genomes
    ## R code:
    pgs_1000g_liab <- pgs_1000g_obs$SCORESUM * sqrt(h2o_to_h2l(K = 0.01))

    r2l = var(pgs_1000g_liab)
    ```

After preparing the input, the BPC function can be applied as follows to compute the predicted disorder probabilities (with a population prevalence of 0.01 and a prior disorder probability of 0.5):

```         
## R code:
source(bpc.R)
pred_prob <- bpc(pgs_liab = pgs_liab, K = 0.01, P = 0.5, r2l = r2l)
```
