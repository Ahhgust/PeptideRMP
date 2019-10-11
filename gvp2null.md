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
For input:
- The predicted protein products from some exome/genome sequencing project. Only [bcftools csq](https://samtools.github.io/bcftools/bcftools.html#csq) (version 1.9-197-g4e51a29) is supported, though more recent versions likely work.
- The protein sequences themselves as provided by ensembl. E.g., Homo_sapiens.GRCh37.pep.all.fa.gz, the v85 file is found here: <br>ftp://ftp.ensembl.org/pub/grch37/release-85/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz
- And a gvp panel file. This file describes the gvps sought. The file must be tab-separated and look like:

| ensembl_protein_id | peptide_seq | start | end | chromosome |
| :----------------: | :---------: | :--:  | :-: | :--------: |
| ENSP00000257198    | AASSQTPTMCTTTVTIK  |  445  |   461  | 18 |
| ENSP00000257198    | AASSQTPTMCTTTVTVK  |      445  |   461 | 18 |
| ENSP00000381072    | ADFSGMSAEK  |      297  |   306  |   18 |
| ENSP00000381072    | ADFSGMSTEK  |     297   |  306   | 18 |

<br>
Notes on the gvp file:
The column order is arbitrary (whatever order you like)<br>
Extra columns are okay! <br>

**But**, make sure you specify the column-names exactly... any typo there will break things.
<br>
Also, pay attention to any warning messages; peptides that were not found will be explicitly presented. The RMP calculator will treat these alleles as being at frequency 5/2n; in fact the whole chromosome will be at this frequency. Several things can cause an allele to fail to be detected in exome sequencing (inaccessible to the array capture is one of them). In general, I recommend you constrain yourself to alleles that have been detected; this can be done with the -i flag to the program. E.g., adding -i AGGSYGFGGAR,ALETLQER would ignore these two alleles (equivalent to deleting them from the panel file)
<br>

Example usage:
> ./gvp2null.py -g gvpPanel.tsv -c AllProteinCodingConsequences.txt.gz -p Homo_sapiens.GRCh37.pep.all.fa.gz -i AGGSYGFGGAR,ALETLQER,GTVEPQLEAR | gzip -9 > protConsequences.tsv.gz


If you:
> zcat protConsequences.tsv | head -n 4

It will look like:

| ID    |  Chromosome  | ProteinID | ProteinStart  |  ProteinStop  |   DetectableAllele1 |  DetectableAllele2 | InferredAllele1 |  InferredAllele2   |
| :---: | :---------:  | :------:  | :-----------: |  :----------: |   :---------------: |  :---------------: | :-------------: |  :------------:    |
| HG00096 | 18  |    ENSP00000257198 | 91 |     101   |  SFSLFLSDGQR  |   SFSLFLSDGQR |    KSFSIFLSDGQR |   KSFSIFLSDGQR |
| HG00100 | 18  |    ENSP00000257198 | 91  |    101   |  SFSLFLSDGQR   |  SFSLFLSDGQR  |   KSFSIFLSDGQR  |  KSFSIFLSDGQR |
| HG00101 | 18  |    ENSP00000257198 | 91  |    101   |  SFSLFLSDGQR   |  nd  |    KSFSIFLSDGQR |   KSFFIFLSDGQR |


<br><br>
Where ID is the Sample ID (from the 1000 genomes project), next is the protein and the protein coordinates, and the two alleles in this individual.
Detected alleles are found by substring search (are invariant to how indels are encoded, which can impact the coordinates of the peptide)
The two alleles are *phased*. I.e., the left-allele is, for example, the maternal, and the right is the paternal (or vice versa).
Sometimes an allele isn't detected; Rather, none of the alleles at a given locus are detected. This becomes "nd" (non-detects).
The allele sequence (including the base-prior) based on the protein coordinates is given as well; this is the inferred allele. Thus we can see that
indivual HG00101 has a nd for their (say) paternal allele; looking at the InferredAllele2 for this individual we can see that there's a KSFF instead of a KSFS prefix to the peptide, which indeed would cause the allele to fail to be detected.

<br>
Special notes:
The detectable alleles are represented in a 19-amino acid alphabet (I -> L), while the original alphabet is shown in the Inferred Allele.
There is also a trypsin bias with the approach; the preceding amino acid of a detected allele must be an R or a K. i.e., the approach assumes a trypsin digest.
Future iterations will make this more flexible...

<br>
As always, if you encounter a bug, contact me!

august DOT woerner AT unthsc DOT edu




