##Prediction and cross-validation

###1. Simple prediction
One can construct a disease risk model based on inference parameters and apply it to separate data using the option 

    -pr

Note that this inference/prediction workflow is meaningful only when the training and test data sets are distinct. For IL, a typical application would be to derive single-SNP parameters for all SNPs in the whole-genome data ("manhattan plot"), and then estimate disease risk of a new set of individuals based on a subset of SNPs. For CL, on the other hand, one would normally derive a subset of SNPs first, infer parameters within this set using the training data, and predict risk probabilities of individuals in the test data using the same SNPs.

####IL

For IL prediction, one reads a result file (option --par file.par) and specifies data files for the test:

    $ gedi -il -pr --par file.par --tped file.tped --tfam file.tfam --pcut 1.0e-7

The parameter file "file.par" is a GeDI output from a separate analysis. Note that parameters must be in IL parameter format. The program reads the parameters and constructs a disease risk model. The option 

    --pcut 1.0e-7

is used to select a subset of SNPs for this process based on a p-value cutoff. The default is 0.05/(no. of SNPs), which corresponds to the Bonferroni correction. Only SNPs with p-values below (or equal to) the cutoff are included in the model. To include all of them, use 

    --pcut 1.0

  The output has the format
	
    0.996    1
    0.892    1
    0.023    0

where each line lists predicted disease risk probability (first column) and phenotype (0=control, 1=case) of the test data. The disease prevalence used in the model is by default read from the parameter file (header line). It can be changed with

    --prev 0.1

*Caution: if the parameters were inferred with 

    -ma_ctl

be sure to set it in prediction as well.

####CL

The CL prediction works similarly:

    $ gedi -cl -pr --par file.par --bfile file --out cl.rsk

which will write the disease risk of each individuals in the test set to "cl.rsk" in the same format as in IL.
However, p-value cutoff is not used and all SNPs in "file.bed" will be included. This is because CL parameter inference result ("file.par") would have been on a pre-selected SNP set, and one would do prediction on a test data ("file.bed" etc) based on the same set of SNPs. If a SNP named in the parameter file cannot be found in the test data file, an error will occur. This search will also assume that the SNPs in parameter file appear in the same order as in the data file. 

Because the SNP selection will be based on parameter file, p-value cutoff is not used and --pcut will be ignored if present. This procedure will also use the alpha parameter read from the parameter file for prediction, which assumes the test set has the same case/control composition as the training set. If a different prevalence is desired for prediction, use the option "--prev 0.1" and the program will recalculate prediction alpha parameter based on the prescribed and training set prevalences. 

###2. Cross-validation

More often, one has a single data set and wishes to perform cross-validation using training and test sub-samples. To do this, simply omit the option "--par file.par":

    $ gedi -il -pr --bfile file --pcut 1e-20
    $ gedi -cl -pr --bfile file --pcut 1e-10 -genotypic

In cross-validation, case and control samples are divided into _k_ ("multiplicity") equal subsets, from which _k_ training/test set pairs are formed:

                case                            control
                |---------------------------|   |------------------------------|

    training	|--------------------|          |-----------------------|
    test                             |------|                           |------|

    training	|-------------|      |------|   |---------------|       |------|
    test                      |------|                          |-------|

    training	|------|      |-------------|   |-------|       |--------------|
    test               |------|                         |-------|

    training	       |--------------------|           |----------------------|
    test        |------|                        |-------|


If the number of individuals _n_ in case/control groups is not a multiple of _k_, the subset size is floor(_n/k_), with the final subset size slightly larger. The multiplicity can be specified using the option

	--cv 5

for which the default is 5. An error is issued if _k_ is too large (subset size less than 2 individuals). The output is a combined list of predicted risk and phenotypes of _k_ test sets (nearly) equal in size to the original data set. Each cross-validation run consists of IL/CL inference on the training set and prediction on the test set. After _k_ runs, GeDI will calculate the area under the curve (AUC) of the receiver operating characteristics (ROC) and print it to standard output along with 95% confidence intervals. Each run will in general have different numbers of SNPs selected under the specified --pcut, because only the test set is used for SNP selection. The standard output will also print the mean number of SNPs selected. In CL, the (geometric) mean of overall asymptotic p-values is also printed. The output file (default "gedi.out") will contain the risk probabilities of the combined test sets. Multiple parameters can be specified for batch runs:

        --ld 0.1 0.2 1
	--eps 1e-3 1e-2 0.1

and cross-validation runs will be repeated. By default, AUC values are output seperately into a file "gedi.auc".

One possibility for the input data is to concatenate tped files from all chromosomes into a single file. Better yet, use a binary ".bed" file. 

***
[Up](README.md)

[Next](tests.md)
