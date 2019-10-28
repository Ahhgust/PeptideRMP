# Files from the FSI:Gen paper:
* EachAllele.tsv
  * Each individual allele encoded so that we can estimate the RMP (sample IDs are unique)
* Pairs.tsv
  * Pairs of alleles that are associated with the same chromosome (heterozygotes are included here; filtered out in the paper)
* sampsToPops.tsv
  * A lookup table associating a sample ID to a population. Used in the Fst computation. **The samples used are an intersection of this file to those found in the consequences prediction.**
  * The samples in the intersection are in: SamplesInIntersection.tsv
* AlleleDetectionsRaw
  * Allele detects in different samples; includes hg38 coordinates (nucleic acid space) and estimates o f the allele frequency from the ExAC project (NFE population).


