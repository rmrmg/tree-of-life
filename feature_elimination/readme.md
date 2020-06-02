# Recursive Feature Elimination with Logistic Regression

For details, see Section S6.2 in Supporting Information. Option `-h` may be helpful as well. Use option `--log name_of_log_file` to write the output to file.

** Note on vectorisation. Order of rdkit descriptors -columns in compressed numpy arrays- is given in `rdkit_desc_names_sorted.txt`.**

## Reproduction of Tables S4 and S5

```
python imb-logistic_rfecv.py --corr_th 0.9 --standarize
python imb-logistic_rfecv.py --corr_th 0.6 --standarize
python imb-logistic_rfecv.py --corr_th 0.5 --standarize
python imb-logistic_rfecv.py --corr_th 0.7 --standarize
python imb-logistic_rfecv.py --corr_th 0.8 --standarize
```

Features obained with `corr_th=0.7` were put in the Table S5. Note that the results may differ slightly due to stochastich character of cross-validation procedure. Names of the features are the names of corresponding functions in `rdkit.Chem.Descriptors` module. Names in the paper where 'translated' according to RdKit documentation.

## Reproduction of Tables S6 and S7

VSA features and framgentts (feature with names starting with `fr_`) are excluded.

```
python imb-logistic_rfecv.py --corr_th 0.9 --standarize --exclude VSA fr_
python imb-logistic_rfecv.py --corr_th 0.6 --standarize --exclude VSA fr_
python imb-logistic_rfecv.py --corr_th 0.7 --standarize --exclude VSA fr_
python imb-logistic_rfecv.py --corr_th 0.8 --standarize --exclude VSA fr_
python imb-logistic_rfecv.py --corr_th 0.5 --standarize --exclude VSA fr_
```
