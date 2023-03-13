# toxin-expression-neotropical-snakes

This repository contains the code used to plot gene expression (FPKM and TPM values) from venom transcriptomes of Neotropical snakes (Colubrinae).

The R-script was run on a local machine where the data (.csv) were stored.

The 'data' directory contains files (.csv) for each venom gland transcriptome containing the sequence ID (output from Trinity), accession number, gene description and e value (output from ncbi BLAST), transcripts per million (TPM) and fragments per kilo million (FPKM). RNA sequences were generated at the University of Michigan Advanced Genomics Core and assembled using Trinity version 2.6.6.

The 'plots' directory contains the resulting figures output by the Rscript. This scripts was developed by Jenna Crowe-Riddell, Kristy Srodawa and Peter Cerda.

The 'outputs' directory contains output data as a table (.csv) that were used to make their respective plots.
