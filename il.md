###Independent loci analysis

Here, a set of genotype/phenotype data are read and each SNP is evaluated separately for association with disease.

    $ gedi -il -tped file.tped --tfam file.tfam -dominant 

will perform genotype distribution inference using the dominant model (DOM). The output has the format

    Chr  SNP   Position  MA   Model   n     alpha   OR   q     p-value   Pd: 0.652463
    1    rs755 194571    T    DOM    3308   0.252  0.604 47.1  6.5e-12
    1    rs107 195861    G    DOM    3308   0.344  0.467 105   7.74e-25

where MA is the minor allele, n is the number of individuals with non-missing genotypes, and q is the likelihood ratio statistic. The column labeled alpha lists contributions each locus makes to the constant term in the exponent of disease risk: alpha = ln[p_d/(1-p_d)]+sum_i(alpha_i) [[Ref.1](wiki/pubs)]. The first line (header) has an extra column recording the disease prevalence p_d. Together, they can be used in a separate analysis to reconstruct alpha. Individual SNP contributions are listed instead of their sum because often one would filter SNPs and include only a part of them to construct disease models. The column "OR" shows the odds ratio contributions of each SNP, which equals exp(beta_i).

If the genotypic model is used (-genotypic), the output has an extra column because there are two ORs: beta_i(a) where a = 1,2.

Logistic regression can be specified for IL. This feature, however, is included for comparison purposes only, and is not recommended: for IL, it performs numerical optimization to obtain solutions that are known analytically from DDA:

    $ gedi -il -lr --tped file.tped --tfam file.tfam

The output format and other possible options are the same as for DDA. Change the tolerance and maximum iterations by 

    --tol 1.0e-8 --itmax 10000

With a sufficiently small tolerance, LR should agree with the analytical DDA result (the default).

***
[Up](README.md)

[Next](cl.md)
