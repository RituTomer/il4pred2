# il4pred2
An updated in silico tool for identification of IL4 Inducing Peptides.

## Introduction
IL4pred2 is a computation based approach to discriminate IL4-inducing peptide in host dependent manner. This tool is able to discriminate between host IL4-inducing peptide from IL4 non-inducing peptides. In the standalone version, we have provided six different models for prediction. These six models are different in terms of their data type and host type. The detailed description of data and host is given below: 
- (i) Host Human, Main Model : This module is developed using MHC class-II binding, IL4-inducing peptides as positive data and both MHC class-II binders and non-binders with IL4 non-inducing property as negative data.
- (ii) Host Human, Alternate 1 Model : This module is developed using MHC class-II binding, IL4-inducing peptides as positive data and MHC class-II binders with IL4 non-inducing property as negative data.
- (iii) Host Human, Alternate 2 Model : This module is developed using MHC class-II binding, IL4-inducing peptides as positive data and MHC class-II non-binders with IL4 non-inducing property as negative data.
- (iv) Host Mouse, Main Model : This module is developed using MHC class-II binding, IL4-inducing peptides as positive data and both MHC class-II binders and non-binders with IL4 non-inducing property as negative data.
- (v) Host Mouse, Alternate 1 Model : This module is developed using MHC class-II binding, IL4-inducing peptides as positive data and MHC class-II binders with IL4 non-inducing property as negative data.
- (vi) Host Mouse, Alternate 2 Model : This module is developed using MHC class-II binding, IL4-inducing peptides as positive data and MHC class-II non-binders with IL4 non-inducing property as negative data.

IL4pred2 is also available as web-server at https://webs.iiitd.edu.in/raghava/il4pred2. Please read/cite the content about the IL4pred2 for complete information including algorithm behind the approach.

## Standalone
The Standalone version of il4pred2 is written in python3 and following libraries are necessary for the successful run:
- scikit-learn
- Pandas
- Numpy

https://webs.iiitd.edu.in/raghava/il4pred2/download.html

## Minimum USAGE
To know about the available option for the stanadlone, type the following command:
```
python il4pred2.py -h
```
To run the example, type the following command:
```
python3 il4pred2.py -i test_sequence.fasta
```
This will predict if the submitted sequences are IL4-inducers or IL4 non-inducers. It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma seperated variables).

## Full Usage
```
usage: 	il4pred2.py [-h] 
                       [-i INPUT ]
                       [-o OUTPUT]
		       [-s host {1, 2}]
		       [-j {1,2,3}]
		       [-t THRESHOLD]
		       [-d {1,2}]
```
```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence(s) in FASTA format
                        or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -s {1,2}, --method {1,2}
                        Host Type: 1: Human, 2: Mouse, by default 1
  -j {1,2,3}, --job {1,2,3}
                        Job Type: 1:Main, 2:Alternate 1, 3:Alternate 2, by default 1
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.5
  -d {1,2}, --display {1,2}
                        Display: 1:Only IL4-inducing Peptides, 2: All
                        peptides, by default 1
```

**Input File:** It allow users to provide input in the FASTA format.

**Output File:** Program will save the results in the CSV format, in case user do not provide output file name, it will be stored in "outfile.csv".

**Job:** User is allowed to choose between different modules, such as, 1 for main model, 2 for alternate 1 model, 3 for alternate 2 model, by default its 1.

**Threshold:** User should provide threshold between 0 and 1, by default its 0.5 for other jobs except PQ Density for which the default threshold is 0.5.

**Display type:** This option allow users to fetch either only il4-inducing peptides by choosing option 1 or prediction against all peptides by choosing option 2.

IL4pred2 Package Files
=======================
It contains following files, brief description of these files given below

INSTALLATION                        : Installations instructions

LICENSE                             : License information

README.md                           : This file provide information about this package

model.zip                           : This zipped file contains the compressed version of model

Feature.zip			    : This zipped file contains the best features used to develop model

pfeature_standalone.zip 	    : This zipped file contains the feature standalone to calculate features

il4pred2.py                          : Main python program


# Reference


