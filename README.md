# Peptide RMP

The code herein can be used to estimate the [random match probability](https://en.wikipedia.org/wiki/Random_match_possibility) (RMP) of set of peptides detected by [liquid chromatography tandem mass spectrometry](https://en.wikipedia.org/wiki/Liquid_chromatography%E2%80%93mass_spectrometry) (LC-MS/MS).
<br><br>
As a refresher, the RMP is used to quantify the rarity of a set of genetic markers. In modern parlance *genetic* is often conflated/equated with/to DNA, but the meaning of genetic instead relates to heredity (this make sense as genetics as a discipline existed long before the discovery of DNA). In fact, some of the earliest studies in (forensic) population genetics did not use DNA markers, but instead were based on allozyme variation. Given this, it is not surprising that peptide markers can be used in forensics, but the *how* is an open question.

## Some methodological considerations

Some problems that we need to consider: 
1. **Peptide markers are not independent**. Typical forensic markers are constructed to be far apart in units of [centimorgans](https://en.wikipedia.org/wiki/Centimorgan), which argues for their (near) statistical independence. This is very convenient; if you compute the RMP of a marker, you can then apply the product rule across markers to estimate the RMP for a set of markers. Peptides detected by LC-MS/MS may come from proteins that are close (in the genome). Alleles that are close (in genetic units) are likely co-inherited, thus the use of the product rule on these alleles will likely overstate the RMP (biasing it downwards).
2. **Peptide markers may be biased towards particular alleles**. Starting with an example, if we only detect a single allele *A*, we don't know if the individual was homozygous *A/A*, or heterozygous *A/B* for some other allele *B*. Further, it may be that the assay is incapable of detecting the *B* allele, or it could also be that the *A* and *B* alleles are differentially expressed. Thus a lack of a detection cannot be considered a random "drop out" event; instead this lack of detection can be the result of any of a slue of stochastic processes that yielded the given result. Trying to interpret any lack of detection should be viewed with a healthy amount of skepticism and for forensic purposes a "drop out" parameter to the RMP is inadvisable.
3. **The rarity of a peptide marker is not the same as the rarity of a particular SNP**. Most commonly the rarity of a peptide is assessed using the allele frequency of some non-synonymous SNP. If the SNP allele is at frequency, say, 10%, then the peptide associated with that allele is assumed to be at the same 10% frequency. This of course is incorrect; other SNPs in the same codon may generate the same peptide (making the peptide frequency >10%), other protein-changing SNPs in the same peptide region may lead to a non-detect (as they change the mass of the peptide), and still other events (e.g., frame-shift mutations) that are perhaps far away may change the reading frame and lead to a lack of peptide detection (again, reducing the allele frequency of some observed peptide). Thus a better estimator of a peptides' frequency is not one based on a single SNP (which projects a peptide from the protein space to the nucleic acid space), but is instead based on the constellation of nucleotide variants found in some individual's haploid genomes. A good way of estimating the rarity of a peptide is to use the predicted proteins from some genomic assay, for example using [bcftools csq](https://samtools.github.io/bcftools/bcftools.html#csq) to predict not just which SNPs are synononymous and non-synonymous, but full protein products that consider the [phase](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007308) of the alleles in question. Using protein predicitions we can project, say, whole exomes into their expected protein products (projecting DNA into protein), as opposed to projecting peptides into DNA (which is, strictly speaking, not possible).
4. **The RMP is not the same as the allele frequency**. In forensic applications allele frequencies are estimated in individuals sampled from a population. Human populations are oftened structured; one type of structure occurs when you sample some "population", but instead what you're measuring is instead two populations that share some alleles. From this its pretty straight-forward to see that an estimated allele frequency (which was estimated not knowing this structure) will not give you a true estimate of the rarity of this allele in the (sub)populations. To address this (and other issues like inbreeding and small deviations from a lack of independence), forensic match statisics use a *&Theta;* correction

## Overview of approach

The basic idea of the RMP approach herein is as follows: 


