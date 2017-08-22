# dgSeq

This repository provides a reference implementation of dgSeq as described in the manuscript:

Disease Gene Prediction by Integrating PPI Networks, Clinical RNA-Seq Data and OMIM Data

Ping Luo, Li-Ping Tian, Jishou Ruan, and Fang-Xiang Wu

dgSeq predicts disease genes from a PPI network and a differential network constructed by RNA-Seq data.

Basic Usage

Example

Input

The input files include the PPI network, the differential network (or the matrix P) and the benchmark disease genes and non-disease genes

Output
Algorithms in the directory "dgSeq/Algorithms/" can compute the AUC values of the prediction as well as the top 100 unknown genes.


Please send any questions you might have about the code and/or the algorithm to ping.luo@usask.ca or faw341@mail.usask.ca
