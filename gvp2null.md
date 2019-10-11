# Making nulls!

This describes how to run the program: gvp2null.py


## Program requirements
This script is written in Python (3.*), and it has no special libraries.
<br>
First, choose a version of Ensembl, and use this version throughout all analyses.
Any and all files provided (canned data) are generated using Ensembl v85, GRCh37
(which is the same as the ExAC project)
<br>

###
For input: <br>
- The predicted protein products from some exome/genome sequencing project. Only [bcftools csq](https://samtools.github.io/bcftools/bcftools.html#csq) (version 1.9-197-g4e51a29) is supported, though more recent versions likely work.
<br>
- The protein sequences themselves as provided by ensembl. E.g.,
<a href="ftp://ftp.ensembl.org/pub/grch37/release-85/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz">Homo_sapiens.GRCh37.pep.all.fa.gz</a>
(v85)
<br>

A gvp panel file. This file describes which peptide alleles are found where...



