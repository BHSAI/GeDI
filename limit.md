##Current limitations

1. IL prediction, including cross-validation, cannot use binary file input:

    -il -pr --bfile file

2. IL cross-validation cannot use logistic regression:

    -il -pr -lr
