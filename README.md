
# Implementation of the Chou-Fasman algorithm  

This is my implementation of the Chou-Fasman method according to ([Chou and Fasman, 1974](https://pubs.acs.org/doi/pdf/10.1021/bi00699a002)).

The structure propensities were adopted from ([Sing et al., 2010](https://pdfs.semanticscholar.org/fd8c/c95aec2d7af19ed28eea3688b3c231d0e745.pdf)) because the original paper 
only evaluated 19 proteins.

## Prerequisites
- Linux or Windows
- Python 3

## Getting started
Clone this repository
```shell
https://github.com/jseidel5/chou-fasman-algorithm.git
cd chou-fasman-algorithm
```
Input
```shell
python algorithm.py sequence_file [propensities_file]
```

sequence_file : Input the name of your *.fasta file in the data directory
containing the sequence.

propensities_file : Input the name of your *.txt file in the data directory 
containing custom propensities.

Note: If no propensities_file is called the default one will be used.

sequence.fasta, helix.fasta and sheet.fasta contain example proteins.