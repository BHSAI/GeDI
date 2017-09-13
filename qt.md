## Quantitative traits

Phenotypes that are quantitative (continuous rather than case-control) can be analyzed by the keyword -qt. The program will expect  file.fam (or file.tfam) file, where the last column (1 or 2 for control or cases) contains the quantitative trait values for each individual. 

### 1. IL

    $ gedi -qt -il --tped file.tped --tfam file.tfam
    $ gedi -qt -il --bfile file
    
The IL mode is similar to that for binary phenotypes. Each SNP read will be analyzed independently and the outcome will be printed (default: gedi.out, can be changed by --out file.out) as a table of odds ratios and p-values. The -qt -il mode uses linear regression for analysis. It differs from case-control data by the model employed: minor allele counts a=0,1,2 for each SNP is used as the predictor in the linear model (a.k.a. additive model: y=beta0+beta*a) by default. If one uses --dominant or --recessive, the predictor is changed into a=0,1.

IL analysis in -qt can include covariates. Covariate files can be read by

    $ gedi -qt -il --bfile file --covar file.covar
    $ gedi -qt -il --bfile file -covar
 
The argument of --covar specifies the covariate file. Without argument (second example), -covar implies that the file is named "file.covar". Error will occur if it doesn't exist. The format of covariate files is the same as in PLINK: it starts with a header line listing the name of covariates after family and individual IDs:

    FID IID COV1 COV2   COV3
    f1  i1  1    0.5    1.1   
    f2  i2  0    0.2    1.0
    f3  i3  0   -0.1    0.8
    
The first two columns must match those of fam (tfam) file. In the above example, there is one discrete and two continuous covariates. These variables are added into the linear model. The GeDI -qt -il output will match those from [PLINK](https://www.cog-genomics.org/plink2) (v1.9) command:

    $ plink-1.90 --bfile file --covar file.covar --linear no-x-sex --hide-covar
  
In other words, sex is not included as a covariate unless it is present in file.covar and male x-chromosome alleles are coded as homozygous.
  
    $ gedi -qt -cl --bfile file --lh 0.1 --ld 0.1 1 10 --pcut 1
  
### 2. CL

In CL, DDA is used with penalizers applied to single-SNP and interaction terms, unless linear regression is explicited requested with -lr as follows:

    $ gedi -qt -cl -lr --bfile file --ld 0.1 1 10
    
In this case, the inference becomes ridge regression (penalized linear regression). After each inference, the variance explained (R = corelation between fitted and actual phenotypes) will be printed. The p-value printed is for the null hypothesis of no association, calculated from the normal distribution applied to Fisher-transformed correlation. Note that these outcomes are all sensitive to penalizer values used. To determine penalizers one has to perform cross-validation (see below). Ridge regression uses single penalizer specified as above with --ld. As in case-control data, the inference will be repeated multiple times, each with the value of penalizer specified (0.1, 1, and 10 in the above example). In addition to the standard output, output file (default is "gedi.out"; can be changed using --out file.out) will list parameter values inferred in the same matrix format as in case-control. The --pcut 1 option dictates all SNPs to be included. Otherwise, the linear regression IL p-value will be used to select SNPs for CL. 
  
Quantitative trait DDA inference will be performed without the -lr flag:

    $ gedi -qt -cl --bfile file --lh 0 0.1 10 10 --ld 0.1 1 10
    
In contrast to case-control, the only inference option is pseudo-likelihood. As in case-control, this inference is parallalized and will scale close to linearly when run in parallel with respect to SNP numbers using MPI-enabled code. Note that two distinct penalizers are used in contrast to ridge regression. These apply to single-SNP and interaction parameters, respectively. The inferene will be repeated for each combination of penalizer values (12 runs in the above).

Cross-validation is performed as in case-control data with the flag -pr:
        
    $ gedi -qt -cl --bfile file --lh 0.01 0.1 1 10 --ld 0.01 0.1 1 10 -pr
    
Under cross-validation, the R^2 becomes the squared correlation between predicted and actual phenotypes for test sets, while inference is performed based on training sets only. The (signed) R is analogous to AUC of case-control inferences and can be maximized with respect to penalizer. In addition to the standard output, a file named "gedi.auc" will be generated with each line listing penalizer values, R, and p-value estimate. The p-values assume that R under the null hypothesis is 0 (not necessarily true). To have a refined estimate of p-value, one can do the following:

    #!/bin/bash
    i=0
    while [ $i -lt 10 ]; do
      $gedi -qt -cl --bfile file --lh 0.01 0.1 1 10 --ld 0.01 0.1 1 10 -pr -null --seed $i
      mv gedi.auc gedi_${i}.auc
      let i=i+1
    done

In the bash script above, we run inferences with phenotype labels permuted to sample the null distribution of R. Having calculated the mean value of R0 = 0.032 from the files gedi_0.auc, gedi_1,auc, etc, 

    $ gedi -qt -cl --bfile file --lh 0.01 0.1 1 10 --ld 0.01 0.1 1 10 -pr --r0 0.032

This inference will print p-values with the null R0 specified.


***
[Up](README.md)

[Next](limit.md)
