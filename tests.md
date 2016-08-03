##Statistical Tests

####1.Single-SNP and interaction tests
In collective inference ([CL](cl.md)), one can additionally perform statistical tests for each locus and pairwise interactions. These tests are based on likelihood ratio statistics derived by comparing the likelihood of the full CL model and that of reduced model where single-SNP parameter of a certain site i or interaction parameters of a given SNP pair ij are pre-set as values corresponding to the case + control pooled sample. Add the option

    -pi

to the CL inference command to trigger the calculation of single-SNP p-values. It will write them to a file named "gedi.pi" one SNP per line. 

These p-values are based on the (infinite sample size) asymptotic chi^2 distribution. Since the null distribution becomes biased under a nonvanishing penalizer value, it is usually necessary to perform null sampling separately, and derive p-values from the statistics of the alternative and null hypotheses. For this purpose, use

    -qi

and GeDI will print the bare statistics values to "gedi.qi." However, for single-SNP tests, these values will often be much larger than their corresponding values under the null hypothesis, and it will be difficult to estimate p-values with a finite set of null distribution sampling data. In this case, the asymptotic p-values can be used as qualitative measures of the relative importance of each SNP under CL inference.

To do the analogous calculation for interactions, use
     
     -qij

and each unique pair of SNPs will be tested, with the results written to "gedi.qij." If you use -pij, asymptotic p-values will be printed to "gedi.pij." In contrast to single-SNP p-values, interaction tests based on null sampling are often feasible. The output format is identical to the parameter output in CL (except the absence of the header line), with single locus (p-values) and pairwise q_ij-statistics as the diagonal and off-diagonal elements, respectively. If you want q_i on the diagonal, use

    -qij -qi

These calculations are parallized; each single-SNP and interaction calculations will be distributed to available cores when run with mpirun.

####2. Null sampling

To estimate interaction p-values, one can first perform the q_ij calculation above, and then sample the null distribution using bootstrap resampling:

     $ gedi -cl --bfile file --pcut 1 -genotypic -boot --seed 1 -qij

where 

    -boot

causes the phenotypes of individuals to be shuffled randomly. The option

    --seed 1

sets the random number seed (positive integer). This process needs to be repeated many times to construct the null distribution of each q_ij value. One can use a shell script varying the argument of --seed like this:

     $!/bin/bash
     n=1
     while [ $n -le 100 ]; do
       $ gedi -cl --bfile file --pcut 1 -boot --seed ${n} -qij
       $ mv gedi.qij gedi_${n}.qij
       let n=n+1
     done

Note the use of "--pcut 1"; this is to insist that GeDI select all supplied SNPs, assumed to be pre-selected and written to "file.bed." Otherwise, GeDI will calculate IL p-values of the SNPs, which will turn out to be close to 1 because of -boot, and reject most of them. 

***
[Up](README.md)

[Next](meta.md)
