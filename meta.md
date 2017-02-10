##Meta-analysis
To increase power, multiple samples can be combined into a single meta-analysis. In discrete discriminant analysis (DDA), the likelihood ratio statistics of each sub-samples are summed into a total, to which statistical tests can be performed. As prerequisites, each sample must have identical sets of SNP lists. Input files can be prepared using PLINK.

In both independen-SNP inference (IL) and collective inference (CL), the command is

    $ gedi -il (-cl) --meta flist.txt
 
with (most of) other options. The file "flist.txt" lists tped/tfam files of each sample sequentially:

    file1.tped file1.tfam
    file2.tped file2.tfam

When using binary files,

    $ gedi -il (-cl) --metab flist.txt

where the "flist.txt" lists

    file1 
    file2

In the binary case, it is assumed that file1.bed, file1.bim, and file1.fam etc exist in the current directory, the order of SNPs as well as allele codings in each sample are the same; errors will occur if not. Please note that PLINK will often recode SNP alleles upon sample changes: e.g., to form sub-groups of individuals listed in files sample1.fam and sample2.fam from a binary file, use

    $ plink-1.90 --bfile file --keep sample1.fam --keep-allele-order --make-bed --out file1
    $ plink-1.90 --bfile file --keep sample2.fam --keep-allele-order --make-bed --out file2
    $ gedi -cl --metab flist.txt 
  
In IL, the output file will list the combined likelihood ratio statistic (q), the degrees of freedom (df), and p-value (P) of samples for which there were sufficient statistics. If a sub-sample leads to an inference failure, that particular sample is skipped and df will be less than the total number of samples times 1 (DOM/REC) or 2 (GEN). If inference fails for all samples, NA results.

In CL, the single-SNP and interaction parameters inferred are averaged over samples weighted by the square-root of effective sample sizes with sample weight = 2/sqrt(1/n_case+1n_control). As in single-sample CL, p-values have to be estimated manually by permutation-resampling. 

***
[Up](README.md)

[Next](limit.md)
