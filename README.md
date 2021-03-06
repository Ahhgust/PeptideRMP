# Peptide RMP

The code herein can be used to estimate the [random match probability](https://en.wikipedia.org/wiki/Random_match_possibility) (RMP) of set of peptides detected by [liquid chromatography tandem mass spectrometry](https://en.wikipedia.org/wiki/Liquid_chromatography%E2%80%93mass_spectrometry) (LC-MS/MS).
## About RMPs

As a refresher, the RMP is used to quantify the rarity of a set of genetic markers. In modern parlance *genetic* is often conflated/equated with/to DNA, but the meaning of genetic instead relates to heredity (this make sense as genetics as a discipline existed long before the discovery of DNA). In fact, some of the earliest studies in (forensic) population genetics did not use DNA markers, but instead were based on allozyme variation. Given this, it is not surprising that peptide markers can be used in forensics, but the *how* is an open question.

## Quick start
-  Use canned data: bcftools csq run on the [phase 3 1000 genomes](https://www.dropbox.com/s/parobd9n91cktv9/AllProteinCodingConsequences.txt.gz?dl=1) (dropbox link)
- Individual -> Population lookup table from [phase 3 1000 genomes](DataForPaper/sampsToPops.tsv). See also AFR.tsv, AMR.tsv, EAS.tssv, EUR.tsv, SAS.tsv for the major super-populations from the 1000 Genomes project as defined [here](https://www.internationalgenome.org/category/population/)
- Download a fasta file of all genes from ensembl: ftp://ftp.ensembl.org/pub/grch37/release-85/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz (copy and paste y'all. Ftp hyperlinks don't appear to work)
- [generate a null distribution of rmps](gvp2null.md)
 *Pay attention to the warning messages! I recommend removing all alleles not detected in the 1KG project*
- Compute RMPs. Note: Two estimators the RMP are of use: the naive_rmp (no theta correction) and the rmp (theta corrected). The latter is presented "raw", and may exceed 1.0.

Example usage:
> /usr/bin/python3 ./pRMP.py  -q sampsToPops.tsv -p EachAllele  -c proteinConsequences.1kg.tsv.gz > EachAllele.RMP

## Dependencies:
- numpy
- standard python libraries: itertools, namedtuple, defaultdict, random, argparse, os.path, sys, csv, gzip


## Some methodological considerations

Some problems that we need to consider: 
1. **Peptide markers are not independent**. Typical forensic markers are constructed to be far apart in units of [centimorgans](https://en.wikipedia.org/wiki/Centimorgan) (cM), which argues for their (near) statistical independence. This is very convenient. If you compute the RMP of a marker, you can then apply the product rule across markers to estimate the RMP for a set of markers. Peptides detected by LC-MS/MS may come from proteins that are close. Alleles that are close (in cM) are likely co-inherited, thus the use of the product rule on these alleles will likely overstate the RMP (biasing it downwards).
2. **Peptide markers may be biased towards particular alleles**. Starting with an example, if we only detect a single allele *A*, we don't know if the individual was homozygous *A/A*, or heterozygous *A/B* for some other allele *B*. Further, it may be that the assay is incapable of detecting the *B* allele, or it could also be that the *A* and *B* alleles are differentially expressed. Thus a lack of a detection cannot be considered a random "drop out" event; instead this lack of detection can be the result of any of a slue of stochastic processes that yielded the given result. Trying to interpret any lack of detection should be viewed with a healthy amount of skepticism and for forensic purposes a "drop out" parameter to the RMP is inadvisable.
3. **The rarity of a peptide marker is not necessarily the same as the rarity of a particular SNP**. Most commonly the rarity of a peptide is assessed using the allele frequency of some non-synonymous SNP. If the SNP allele is at frequency, say, 10%, then the peptide associated with that allele is assumed to be at the same 10% frequency. This of course is incorrect; other SNPs in the same codon may generate the same peptide (making the peptide frequency >10%), other protein-changing SNPs in the same peptide region may lead to a non-detect (as they change the mass of the peptide, which lowers the peptide frequency), and still other events (e.g., frame-shift mutations) that are perhaps far away may change the reading frame and lead to a lack of peptide detection (again, reducing the allele frequency of some observed peptide). Thus a better estimator of a peptide's frequency is not one based on a single SNP (which projects a peptide from the protein space to a single allele in the nucleic acid space), but is instead based on the constellation of nucleotide variants found in some individual's haploid genomes. A good way of estimating the rarity of a peptide is to use the predicted proteins from some genomic assay, for example using [bcftools csq](https://samtools.github.io/bcftools/bcftools.html#csq) to predict not just which SNPs are synononymous or non-synonymous, but full protein products that consider the [phase](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007308) of the alleles in question. Using protein predicitions we can project, say, whole exomes into their expected protein products (projecting DNA into protein), as opposed to projecting peptides into DNA (which is, strictly speaking, not possible).
4. **The RMP is not the same as the allele frequency**. In forensic applications allele frequencies are estimated in individuals sampled from a population. Human populations are oftened structured; one type of structure occurs when you sample some "population", but instead what you're measuring is instead two populations that share some alleles. From this it's pretty straight-forward to see that an estimated allele frequency (which was estimated not knowing this structure) will not give you a true estimate of the rarity of this allele in the (sub)populations. To address this (and other issues like inbreeding and small deviations from a lack of independence), forensic match statisics use a *&Theta;* correction, often estimated as F<sub>ST</sub> (see the recent article by [Buckleton et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4899164/) for a description of this). *&Theta;* corrections tend to increase the RMP as well (though not always), which argues that naive approaches to the RMP problem that neglect *&Theta;* may similarly be biased downwards.


## Overview of approach

The basic idea of the RMP approach herein is as follows: 
<br><br>
First, assume that the peptides detected are:
* autosomal (no sex chromosomes)
* single-copy (a given peptide can only stem from a single locus in the genome)
* mostly correct (a low rate of false discovery)
<br><br>
Further, alleles between chromosomes are assumed to be independent and alleles within the same chromosome are permitted any level of dependence. Given this, the code will:
1. **Partition alleles by chromosome.** After which the RMP is assessed at the level of the chromosome (that is, within chromosomes), similar to approaches used with Y chromosome and mitochondrial markers. Products of the RMP (within chromosomes) is taken to yield the final genomic RMP.
2. **Find consistent diplotypes.** Within each chromosome, all pairs of haploid chromosomes that may have yielded the set of observed peptides are found. These are termed the consistent diplotypes, which are constrained to the loci associated with the markers detected. Note that any amount of drop-out is permitted, so long as there is no peptide observed that did not originate from the diplotype.
3. **Estimate the RMPs.** RMPs are estimated based on the consistent diplotypes. As there are many such diplotypes (and many RMPs), the maximum RMP is chosen as the final estimator to ensure that the method is conservative. RMPs are assessed using a *&Theta;* correction  (computed *de novo* as F<sub>ST</sub>) and considering a minumum allele frequency of 5/2*n* (*n* is the diploid sample size).
4. **Simulate drop-in** Drop-in alleles may occur. To model this we use a simple Monte Carlo procedure. Alleles are randomly dropped out of the observed peptides (given some fixed probability *f* applied to all peptides) and the RMP is estimated on the resultant peptides. 

# In practice

The chromosomes used for comparison in practice come from the [1000 genomes project](https://www.internationalgenome.org/category/vcf/). The haploid chromosomes characterized by the project are then converted into proteins using [bcftools csq](https://samtools.github.io/bcftools/bcftools.html#csq) and all pairs of haploid chromosomes are used in the assessment of the RMP. F<sub>ST</sub> is estimated within European populations as found in the 1000 genomes project. 

# Funding

This research is based upon work supported in part by the Office of the Director of National Intelligence (ODNI), Intelligence Advanced Research Projects Activity (IARPA), via contract number 2018-18041000003. The views and conclusions contained herein are those of the authors and should not be interpreted as necessarily representing the official policies, either expressed or implied, of ODNI, IARPA, or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for governmental purposes notwithstanding any copyright annotation therein.


