# BPC

The Bayesian polygenic score Probability Conversion (BPC) approach transforms Bayesian polygenic scores (PGSs) into probabilities of disease. The BPC approach is described in the following paper: [TBA].

### Getting Started

The following files and tools will need to be downloaded:

-   The `bpc.R` and `h2o_to_h2l.R` files. Alternatively you can clone this whole repository using the following git command: `git clone https://github.com/euffelmann/bpc.git`

-   1000 Genomes population reference files (from e.g. <https://ctg.cncr.nl/software/magma>)

-   PRScs (<https://github.com/getian107/PRScs/tree/master>) or GCTB for SBayesR (<https://cnsgenomics.com/software/gctb/#Download>)

### Applying the BPC function

The BPC function takes four inputs:

-   **pgs_liab**: A single individual's PGS (or a vector of PGSs) on the liability scale based on the posterior mean betas from a Bayesian PGS method. While in theory any Bayesian PGS method that is well-calibrated for continuous traits can be used, we have only evaluated PRScs and SBayesR.

    -   PRScs: Use the effective sample size (neff) of the GWAS training sample as input. Using a PGS for schizophrenia as an example (population prevalence = 0.01), transform the PGS from the observed scale with 50% case ascertainment (because neff was used as input in PRScs) to the liability scale with the following code:

        ```         
        source(h2o_to_h2l.R)
        pgs_liab <- pgs_obs * sqrt(h2o_to_h2l(K = 0.01, P = 0.5, h2o = 1))
        ```

    -   SBayesR: To achieve well-calibrated PGSs, the GWAS betas and standard errors need to be transformed to the standardized observed scale with 50% case ascertainment first before SBayesR can be used.

        ```         
        gwas$b5050 <- gwas$b / (se * sqrt(neff))
        gwas$se <- 1 / sqrt(neff)
        ```

        After the GWAS betas and standard errors have been transformed in this way, you can use SBayesR with neff as input and rescale the PGS to the liability scale in the same way as for PRScs.

-   **K**: The population prevalence of the disorder of interest (e.g. 0.01 for schizophrenia)

-   **P**: The prior disorder probability. This can be interpreted as the proportion of cases in a (hypothetical) testing sample.

-   **r2l**: The explained variance (*R*<sup>2</sup>) on the liability scale. In the BPC approach, the variance of a well-calibrated PGS in a population reference sample (e.g. 1000 Genomes) is used to estimate *R*<sup>2</sup><sub>liability</sub>. Simply use the posterior mean betas from PRScs or SBayesR (having used neff as input) to compute PGSs in 1000 Genomes, then transform the PGSs from the observed scale to the liability scale, and compute the variance of the PGSs to estimate *R*<sup>2</sup><sub>liability</sub>. Example code with schizophrenia as the disorder of interest:

    ``` 
    # transform PGSs from observed scale to liability scale in 1000 Genomes
    pgs_1000g_liab <- pgs_1000g_obs * sqrt(h2o_to_h2l(K = 0.01, P = 0.5, h2o = 1))

    r2l = var(pgs_1000g_liab)
    ```

After preparing the input, the BPC function can be applied as follows to compute the probability of disease:

```
source(bpc.R)
pred_prob <- bpc(pgs_liab = pgs_liab, K = 0.01, P = 0.5, r2l = r2l)
```

