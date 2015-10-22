##Collective Inference (CL)

In CL analysis, all single-locus and pairwise interaction effects are inferred simultaneously:

    $ gedi -cl --tped file.tped --tfam file.tfam --pcut 1.0e-7

The input genotype data could either be an intermediate-sized SNP set obtained with a suitable variable selection procedure, or whole genome-data. Because CL analysis is practical only for sets of the order of ~100 SNPs, filtering based on IL p-value cutoff is applied by default. The cutoff value default is 0.05/(no. of SNPs) (Bonferroni correction), which can be modified with 

    --pcut 1.0e-7

Use "--pcut 1.0" to skip filtering (e.g., when using a pre-selected data set).

The DDA will be run for control, case, and then combined ("pooled") groups sequentially. The collective likelihood ratio statistic and (asymptotic chi^2) p-value will be printed to standard output. Resulting parameters are written to the output file in the following format:

    alpha: 0.403 Pd: 0.632
    rs#1  beta1    gamma12  gamma13  gamma14  ... gamma_1m
    rs#2  gamma21  beta2    gamma23  gamma24  ... gamma_2m
    rs#3  gamma31  gamma32  beta3    gamma34  ... gamma_3m
    ...
    rs#m  gamma_m1 gamma_m2 ...                   beta_m

where the first line shows alpha and disease prevalence values. Below the line, the first column lists rs-ID of SNPs, and the rest of columns show the square matrix with diagonal terms given by the single-locus risk parameter beta_i, and off-diagonal terms the pairwise risk parameter gamma_ij. If genotypic model was used, each site has two beta_i, and each site pair has 4 gamma_ij values. These are printed in (2m)x(2m) matrix form (m is no. of SNPs).

There are three methods available for CL:

###1. Pseudo-likelihood (PL)

The default is pseudo-likelihood maximization, specified by

    -pseudo

which can be applied to fairly large data sets without difficulty (the number of SNPs _m_ can be up to hundreds or more). It employs penalizer (lambda1, lambda2), each acting on field h (single-locus) and coupling J (interaction) parameters. They are specified by options such as

    --lh 0
    --ld 0.1

for lambda1, lambda2, respectively (these are defaults). For the latter, multiple values can be specified

    --ld 1e-3 0.01 0.1

to repeat the analysis with different conditions. When multiple values are used, the parameter sets for these values will be concatenated into the output file. This feature is useful when examining the performance of the model with a range of different penalizer values, and can be combined with [cross-validation](CV).

Tolerance and maximum iteration are changed by 

    --tol 1.0e-5
    --imax 10000

(these are defaults), respectively. Use these when the inference appears to hang or numerical problems are suspected.

The PL option is parallelized; when run on multiple processors with mpirun, the scaling will be close to linear.

###2. Mean field (MF)

The mean field option is specified by

    -mf

This option can be applied to large data sets (_m_ up to thousands) with some compromises in bias and variance. It employs a regularization parameter epsilon = [0,1] that controls the relative importance of interactions. The IL limit is reached with epsilon = 0, and epsilon = 1 corresponds to the full interaction limit. It is specified with 

    --eps 0.1 

(default). Although CL should become IL with "--eps 0", the outcomes will not be exactly same because in CL, missing genotypes are counted as major alleles, and MF includes one prior count of individual with uniform genotype frequencies to avoid singularity. One can experiment with different values or use [cross-validation](CV) to find an optimal condition.

###3. Exact enumeration (EE)

One may choose to do an exact enumeration (EE) DDA by specifying the option 

     -ee 

This is the most accurate method but the computational cost scales exponentially with the number of SNPs _m_; execution time will become impractical for _m_ much larger than ~20. Use a p-value cutoff (e.g., "--pcut 1e-7") to reduce the SNP number. EE is useful when trying to validate other approximate treatments (PL or MF). It uses the same penalizing parameters (lambda1, lambda2) as PL with the same syntax. Site and interaction [tests](Tests)
can also be performed to calculate partial p-values of loci and interactions. These values from EE, if it is feasible to calculate them, will be more accurate than PL.

###4. Logistic regression (LR)

One can use logistic regression for CL inference by specifying 

    -lr

This option behaves the same way as PL/EE (also using lambda). LR will in general be faster than DDA but less accurate.

***
[Up](readme.md)

[Next](cv.md)
