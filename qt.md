## Quantitative traits

(Under development)

Phenotypes that are quantitative (continuous rather than case-control) can be analyzed by the keyword -qt

    $ gedi -qt -il --tped file.tped --tfam file.tfam
    $ gedi -qt -cl --bfile file --ld 0.1 1 10 --pcut 1
  
The program will read file.fam (or file.tfam) file, where the last column (1 or 2 for control or cases) contains the quantitative trait values for each individual. In IL, linear regression is used to each SNP independently. In CL, DDA is used with penalizer applied to single-SNP and interaction terms, unless linear regression is explicited requested with -lr as follows:

    $ gedi -qt -cl -lr --bfile file --ld 0.1 1 10
    
In this case, the inference becomes ridge regression (penalized linear regression). After each inference, the variance explained (R2 = corelation between fitted and actual phenotypes squared) will be printed. The --pcut 1 option dictates all SNPs to be included. Otherwise, the linear regression IL p-value will be used to select SNPs for CL.

Cross-validation is performed as in case-control data with the flag -pr:
        
    $ gedi -qt -cl --bfile file --ld 0.1 1 10 -pr --pcut 1
    
Under cross-validation, the R2 becomes the squared correlation between predicted and actual phenotypes for test sets, while inference is performed based on training sets only. The R2 is therefore analogous to AUC of case-control inferences and can be maximized with respect to penalizer.
    
    

***
[Up](README.md)

[Next](limit.md)
