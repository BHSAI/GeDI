#!/bin/bash
data='chr19sim'

# IL analysis using tped/tfam file
../../gedi -il --tped ../${data}.tped --tfam ../${data}.tfam --out il_dda1.dat

# IL analysis under genotypic model
../../gedi -il --tped ../${data}.tped --tfam ../${data}.tfam -genotypic --out il_dda2.dat

# IL analysis using logistic regression
../../gedi -il -lr --tped ../${data}.tped --tfam ../${data}.tfam --out il_dda3.dat

# IL prediction using parameter file (training set == test set; don't do this in real work)
../../gedi -il -pr --tped ../${data}.tped --tfam ../${data}.tfam --par il_dda1.dat --pcut 0.1 --out il_pr.dat

# IL prediction with cross-validation
../../gedi -il -pr --tped ../${data}.tped --tfam ../${data}.tfam --pcut 0.1 --out il_pr_cv.dat

# CL analysis with MFA
../../gedi -cl --tped ../${data}.tped --tfam ../${data}.tfam --pcut 0.1 --out cl_dda1.dat 

# CL analysis with EE
../../gedi -cl -ee --tped ../${data}.tped --tfam ../${data}.tfam --pcut 0.1 --out cl_dda2.dat 

# CL with LR
../../gedi -cl -lr --tped ../${data}.tped --tfam ../${data}.tfam --pcut 0.1 --out cl_lr.dat 

# CL prediction (again training==test; no good)
../../gedi -cl -pr --par cl_dda1.dat --tped ../${data}.tped --tfam ../${data}.tfam --out cl_pr1.dat

# CL prediction with different prevalence
../../gedi -cl -pr --par cl_dda1.dat --tped ../${data}.tped --tfam ../${data}.tfam --prev 0.2 --out cl_pr2.dat

# CL cross-validation with MFA
../../gedi -cl -pr --tped ../${data}.tped --tfam ../${data}.tfam --pcut 0.1 --out cl_pr3.dat
