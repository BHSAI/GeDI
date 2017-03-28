##Statistical Tests

####1.Meta-analysis
To increase power, multiple samples can be combined into a single meta-analysis. In discrete discriminant analysis (DDA), the likelihood ratio statistics of each sub-samples are summed into a total, to which statistical tests are performed. As prerequisites, each sample must have identical sets of SNP lists. Input files can be prepared using PLINK.

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

In the binary case, it is assumed that file1.bed file1.bim file1.fam etc exist in the current directory and error wiill occur if not.


***
[Up](README.md)

[Next](limit.md)
