MIsoCor
==============




Introduction
-------
MIsoCor is a software tool for the correction for natural abundance and isotopic impurity in isotope labeling experiments. It is revised from IsoCor, which was published a few years ago. The purpose of the revision is

1. To fix the issues that IsoCor has with correction matrix construction and constrained linear regression. 
2. To provide a MATLAB tool that can easily be incorporated into exisiting pipelines for metabolic tracer and flux analysis. 




Tutorial
-------
There are five MATLAB files in the repository. The MATLAB function file "misocor.p" can be called to perform the correction aforementioned. The example file "misocor_example.m" illustrates how to call the function using experimentally measured fractional abundances of lactate. The problem of IsoCor with the construction of the correction matrix is shown in "isocor_purity_error.m". Using MATLAB symbolic toolbox, it compares the IsoCor approach to construct the a single matrix for both natural abundance and isotopic impurity and the standard approach to construct matrices for natural abundance and isotopic impurity respectively. Since the built-in function of convolution is only numerical, a symbolic version of it, "symconv.m' is attached. The function "parse_formula.m" is written by Jeffrey Kantor, and parses the chemical formula into a struct. 

Though there exist fundamental issues with IsoCor, its contribution to the field is indelible. When there are more than two stable isotopes, explicit expression of the correction matrix is difficult to obtain through combination probabilities because all possible isotopes of an element need be considered in a single matrix. The probelm is solved in an elegant manner in IsoCor using iterative convolution. For anyone who uses this tool, please cite the original IsoCor application note as well. 




Citations
--------
1. Millard, P., et al. IsoCor: correcting MS data in isotope labeling experiments. Bioinformatics 2012;28(9):1294-1296
