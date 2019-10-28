#!/usr/bin/python3
# Written by August Woerner

# MIT License

# Copyright (c) [2019] [August E. Woerner]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



from math import log, exp, factorial
import numpy
from collections import namedtuple
from collections import defaultdict
import sys
import csv
import gzip
from itertools import combinations
import random
import argparse
import os.path


# poor man's #define
DEFAULT_SIM_SAMPLE_SIZE = 50
DEFAULT_REP_SAMPLE_SIZE = 1000
ALT_MIN_THETA = 0.00
NUCLEIC_ACID = 1
PROTEIN = 0

#GenomeCoordinate = namedtuple("GenomeCoordinate", ["chromosome", "position"])
Allele = namedtuple("Allele", ["basecall", "frequency"])

# holds genomic coordinates that correspond to either nucleic acids (chromosoem and position are defined)
# or proteins (chromosome, protStart, protStop and proteinID are defined)
class GenomeCoordinate:
  def __init__(self, chromosome, position=None, protStart=None, protStop=None, proteinID=None):
    self.chromosome = chromosome
    if (position is not None):
      self.position = int(position)
      self.isProt=False
    else:
      self.protStart= int(protStart)
      self.protStop = int(protStop)
      self.proteinID = proteinID
      self.isProt=True
    
  def __hash__(self):
    if self.isProt:
      return hash( (self.chromosome, self.protStart, self.protStop, self.proteinID) )      
    else:
      return hash( (self.chromosome, self.position) )
      
  def __eq__(self, other):
    if self.isProt:
      if other.isProt:
        return self.chromosome==other.chromosome and \
          self.proteinID == other.proteinID and \
          self.protStart == other.protStart and \
          self.protStop == other.protStop
      return False
    elif not other.isProt:
      return self.chromosome ==other.chromosome and \
        self.position == other.position
    return False
  
  def __ne__(self,other):
    return not(self==other)

  def __getitem__(self, key):
    if key == 0:
      return self.chromosome
    if key == 1:
      if self.isProt:
        return self.protStart
      return self.position
    
    return None
  
  def __lt__(self, other):
    if self.chromosome != other.chromosome:
      return self.chromosome < other.chromosome
    # ensure the types are comparable
    if self.isProt + self.isProt == 1:
      return 1
    
    if self.isProt:
      if self.proteinID != other.proteinID:
        return self.proteinID < other.proteinID
      return self.protStart < other.protStart
    # DNA case
    return self.position < other.position
  
  def isProt(self):
    return self.isProt
  def isNuc(self):
    return not(self.isProt)

  def __str__(self):
    if self.isProt:
      return "chrom: " + self.chromosome + "\tProt: " + self.proteinID + "\tCoords: " + \
        str(self.protStart) + "\t" + str(self.protStop)
    return "chrom: " + self.chromosome + "\tCoord: " + str(self.position)
  
  def __repr__(self):
    return(self.__str__())
  
def rmpHomozygousWithTheta(theta, pu):
  """
  Notation from Baldings and Nichols (1994)
  taken from Buckleton FSI Gen (23) 2016; 91-100
  https://doi.org/10.1016/j.fsigen.2016.03.004
  """
  # 2theta + (1-theta)*pu
  numerator1 = (theta+theta) + (1-theta)*pu
  # 3theta + (1-theta)*pu == theta + numerator1
  numerator2 = theta + numerator1
  
  # 1 + theta TIMES 1 + 2theta
  denom = (1+theta) * (1 + theta + theta)

  return (numerator1*numerator2)/denom

def rmpHeterozygousWithTheta(theta, pu, pv):
  """
  Notation from Baldings and Nichols (1994)
  taken from Buckleton FSI Gen (23) 2016; 91-100
  https://doi.org/10.1016/j.fsigen.2016.03.004
  """

  # theta + (1-theta)pu
  numerator1 = (theta + (1-theta)*pu)
  numerator2 = (theta + (1-theta)*pv)

  denom = (1+theta)*(1+theta+theta)

  return (2*numerator1*numerator2)/denom


# pure python n-choose-k implementation
def nChooseK(n, k):
  return factorial(n)//factorial(k)//factorial(n-k)


def isoleucine2leucine(s):
  """
  converts all isoleucines to leucines in a string...
  this is done to remove equivalent weight amino acids
  it serves to reduce the alphabet from 20 aa to 19.
  """

  ss = list(s)

  for i in range(len(ss)):
    if ss[i] == 'I':
      ss[i] = 'L'

  return "".join(ss)


def organizeProt(proteomeHitsFile, mode=NUCLEIC_ACID, protLUTFile=None, autosomesOnly=True, missenseOnly=False):
  """
  This function parses a *rmp* file from sigsci
  and creates a datastructure that associates:
  id -> genome coordinates -> allele(s)
  The outer two datastructures are dictionaries, and the inner key (alleles) is a list
  """

  expectedPeps = {}
  if mode != NUCLEIC_ACID:
    with open(protLUTFile, "r") as fh:
      reader = csv.DictReader(fh, delimiter="\t")
      for row in reader:
        prot = row["ProteinID"]
        pep = isoleucine2leucine(row["Peptide"])
        if pep not in expectedPeps:
          expectedPeps[ pep ] = (prot, row["Start"], row["Stop"]) # start and stop coordinates as defined in the expected peptides file
        elif expectedPeps[pep] != (prot, row["Start"], row["Stop"]):
          print("Inconsistent coordinates detected for peptide: " , pep , \
                (prot, row["Start"], row["Stop"]) , \
                expectedPeps[pep], \
                file=sys.stderr)
          sys.exit()
          
  prots = {}
  
  with open(proteomeHitsFile) as tsv:
    reader = csv.DictReader(tsv, delimiter="\t")
    
    for row in reader:
      if autosomesOnly:
        if row["chromosome"] == 'chrX':
          continue
        if row["chromosome"] == 'chrY':
          continue
        if row["chromosome"].startswith('chrM'): # can be chrM or chrMT (thank you ensembl)
          continue
      if missenseOnly and "allele_change_type" in row and  row["allele_change_type"] != 'missense':
        continue

      if "label" in row:
        ID = row["label"].split(" ")[-1] # the last element is the sample name. kinda ugly...
      elif "donor_id" in row:
        ID = row["donor_id"]
      elif "sample_id" in row:
        ID = row["sample_id"]
      else:
        print("Can't find the column for the sample ID!", file=sys.stderr)
        sys.exit(1)
        
      # added by AW.
      # this adds support for peptide allele calling (without hand-waving!)
      if mode == PROTEIN:
        peptide = isoleucine2leucine(row["peptide_seq"])
        if peptide not in expectedPeps:
          print("Skipping allele: ", row["peptide_seq"] ,\
                " it wasn't typed and/or it was skipped in the genomic dataset (pep not found)", \
                file=sys.stderr, sep="")
          continue
        
        protStart = int( expectedPeps[peptide][1])
        protEnd = int( expectedPeps[peptide][2])
        prot = expectedPeps[peptide][0]

        chrom = row["chromosome"] # TODO: move this to the LUT!
        f = chrom.find(";") # sometimes it's a ; separated list. not fun!
        if f != -1:
          chrom = chrom[:f]
          
        if chrom.startswith("0"):
          chrom = chrom[1:]
        elif chrom.startswith("chr"): # force ensembl chromosome names. 
          chrom = chrom[3:]
        # python doesn't support overloaded constructors. boo. hence the none.
        coord = GenomeCoordinate(chrom, None, protStart, protEnd, prot)
        # TODO: generalize!
        allele = Allele(isoleucine2leucine(peptide), -1)
        if ID not in prots:
          prots[ID] = {}
          
        if coord not in prots[ID]:
          prots[ID][coord] = [allele]
        else: # add support for heterozygous calls
          gotIt=False
          for a in prots[ID][coord]:
            if a[0] == allele[0]: # duplicate entry. Not nice, but okay
              gotIt=True
              break
              
          if not gotIt:
            alleles = prots[ID][coord].append(allele)
        
        continue

      # handinlg alleles in the nucleic acid space is more annoying. it's a 1:many relationship (peptide -> snps)
      # if multiple SNPs are spanned then the data are ; separated.
      # was tokenized off of gene_allele. That's the old nomenclature
      geneAlleles = row["nucleotide_observed"].split(";")
      
      if "allele_freq" in row:
        freqs = row["allele_freq"].split(";")
      else:
        freqs = [-1] * len(geneAlleles)


      if "nuc_change" in row: # encoded as fromallele_chromosome_hg38start_toallele
        nucChanges = row["nuc_change"].split(";") # may encode multiple variants...
      else:
        if "hg38stop" not in row or "chromosome" not in row:
          print("I need either a nuc_change column or an hg38start and chromosome column!", file=sys.stderr)
          sys.exit(1)
          
        nucChanges = ["", row["chromosome"], str(int(row["hg38stop"])-1), ""]
        if nucChanges[1].startswith("0"):
          nucChanges[1] = "chr" + nucChanges[1][1:]
        elif not nucChanges[1].startswith('chr'):
          nucChanges[1] = "chr" + nucChanges[1]

        nucChanges = [ "_".join(nucChanges) ]

      # each peptide is associated with 1+ genomic alleles
      for nucChange in nucChanges:
        # frombase, chromosome, position, tobase
        nuc = nucChange.split("_")
        coord = GenomeCoordinate( nuc[1], int(nuc[2]))
        allele = Allele( geneAlleles.pop(0), float(freqs.pop(0)) )
      
        if ID not in prots:
          prots[ID] = {}

        if coord not in prots[ID]:
          prots[ID][coord] = [allele]
        else:
          gotIt=False
          for a in prots[ID][coord]:
            if a[0] == allele[0]: # duplicate entry. Not nice, but okay
              gotIt=True
              break
              
          if not gotIt:
            alleles = prots[ID][coord].append(allele)

  return prots

def marginalRMP(a1, a2, theta=0.03):
  """
  This takes the diploid genotype call from the proteome 
  as two Allele namedtuples
  and computes a RMP of just this single genotype call
  This DOES address population structure (theta)
  but theta needs to be estimatated a priori
  """

  if (a1.basecall == a2.basecall): # inferred homozygote
    freq = a1.frequency
    rmpHomo = rmpHomozygousWithTheta(theta, freq)
    rmpHetero = rmpHeterozygousWithTheta(theta, freq, 1-freq)
    return rmpHomo + rmpHetero
  
  freq1 = a1.frequency
  freq2 = a2.frequency
  return rmpHeterozygousWithTheta(theta, freq1, freq2)



def marginalRMPNaive(a1, a2):
  """
  This takes the diploid genotype call from the proteome 
  as two Allele namedtuples
  and computes a RMP of just this single genotype call
  This does NOT address population structure.
  """
  if (a1.basecall == a2.basecall): # inferred homozygote
    freq = a1.frequency

    rmpHomozygous = freq*freq #is a homozygous sample
    rmpHeterozygous = 2 * freq * (1-freq) # is a heterozygous sample (exome) and we're just detecting one allele
    # mutually exclusive possibilitiyies (het, hom), either can contribute to measured phenotype
    # thus you sum...
    return rmpHeterozygous + rmpHomozygous
  else:
    freq1 = a1.frequency
    freq2 = a2.frequency
    return 2 * freq1 * freq2


def getMinRmpByChrom(rmps):
  """
  This takes the output of: getAllMarginalRMPs
  and computes the min RMP for each marker organized by
  chromosome
  The user may use this, taking the product of the RMPs
  however this argues for biological independence of the chromosomes
  which may or may not hold for protein sequences.
  """

  prevChrom = ""
  bestRmps = {} # associates a chromosome with the best RMP
  for coordinate, rmp in rmps.items():

    if coordinate.chromosome not in bestRmps:
      bestRmps[coordinate.chromosome] = (coordinate, rmp)
    elif bestRmps[coordinate.chromosome][1] > rmp: # want the smallest RMP by chromosome
      bestRmps[coordinate.chromosome] = (coordinate, rmp)

  d = {}
  for c,r in bestRmps.values():
    d[c] = r

  return d



def getAllMarginalRMPs(d, lut=None):
  """
  This takes in a datastructure parsed by organizeProt
  (the inner dictionary of it)
  and it produces the RMP of every marker in it.
  this is the RMP of just this one marker
  it cannot be combined by products
  unless some steps are taken to make it independent
  """
  out = {}
  for coord, alleles in d.items():
    if lut is None or coord in lut:
      if len(alleles) == 1:
        rmp = marginalRMP(alleles[0], alleles[0])
      else:
        rmp = marginalRMP(alleles[0], alleles[1])
      out[coord] = rmp
      
  return out


def organizeSnpToSap(snp2sapFile):
  """
  This parses the snp2sap database, and creates a dictionary of lists
  in it a coordinate is associated with the Alleles (Allele namedtuples) at that coordinate
  the 0th Allele in this list is the reference sequence
  Because indel polymorphisms are... annoying, it is possible for the reference allele
  to be repeated as a variant allele (index >0). This happens when insertions and deletions overlap.
  To be safe, the allele forms should be ignored until something more sophisticated is done
  to work with indel polymorphisms.
  In short, the allele frequencies in this case are correct, but the variants may not be.
  """
  if snp2sapFile.endswith(".gz"):
    fh = gzip.open(snp2sapFile, 'rt')
  else:
    fh = open(snp2sapFile, "r")

  
  reader = csv.DictReader(fh, delimiter="\t")

  snp2sap = {}
  
  for row in reader:
    coord = GenomeCoordinate( "chr" + row["chromosome"], int(row["hg38start"]) )

    a2 = Allele( row["allele_polymorphism"], float(row["allele_freq_alt"]))
    if coord not in snp2sap:
      a1 = Allele( row["allele_reference"], float(row["allele_freq_ref"]))
      snp2sap[coord] = [a1, a2]
    else: # multiple variants
      snp2sap[coord].append(a2)
      
  fh.close()


  for coord in snp2sap.keys():
    af = 0
    for a in snp2sap[coord]:
      af +=  a.frequency

    # imprecision can cause the frequency > 1
    # < 1 is okay (there's another allele here < 1% frequency in Europeans...)
    if af > 1:
      alleles = snp2sap[coord]
      newAlleles = []
      # adjust the allele frequencies so that they sum to something <= 1
      for a in alleles:
        newAlleles.append( Allele(a[0], a[1]/af) )
        
      snp2sap[coord] = newAlleles
  
  return snp2sap



def makeHaploidDiploid(d):
  """
  Takes a dictionary of calls that are haploid key/value pairs
  (individual -> haploid genotype)
  and creates all pairs of diploid genotypes
  A smarter construction would use a generator...
  """
  individuals = sorted(d.keys()) # need to ensure that the individuals are deterministically pulled from the dict
  
  out = {}
  for i1,i2 in combinations(individuals, 2): # all pairs
    out[i1 + "\t" + i2] = (d[i1], d[i2])

  return out


def loadGenotypes(filename, siteList, samps2pops, haploid=True, space=NUCLEIC_ACID):
  """
  This parses a VCF/bcftools csq file (converted to tabular form)
  and extracts the genotypes from it that correspond to 
  SNPS found in siteList
  a dictionary of dictionaries is returned.
  the outer dictionary is keyed on the coordinate
  and the inner dictionary associates a sample ID (str)
  to their diploid genotypes (tuple)
  """

  if filename.endswith('gz'):
    fh = gzip.open(filename, "rt")
  else:
    fh = open(filename, 'r')
    
  reader = csv.DictReader(fh, delimiter="\t")

  dictOfDicts = {}

  if space == NUCLEIC_ACID: 
    for row in reader:

      if row["ID"] not in samps2pops:
        continue
      
      if row["Chromosome"].startswith('chr'):
        coord = GenomeCoordinate(row["Chromosome"],  int(row["hg38start"]) )
      else:
        coord = GenomeCoordinate("chr" + row["Chromosome"],  int(row["hg38start"]) )
    
    
      if coord in siteList:
        if coord not in dictOfDicts:
          genotypes = {}
          dictOfDicts[coord] = genotypes
        else:
          genotypes = dictOfDicts[coord]
        
        if haploid:
          genotypes[row["ID"]+"_A"] = row["Allele1"]
          genotypes[row["ID"]+"_B"] = row["Allele2"]
        else:
          genotypes[row["ID"]] = (row["Allele1"], row["Allele2"])
  else:
    for row in reader:

      if row["ID"] not in samps2pops:
        continue
      
      # some chromosomes are encoded with leading 0s... no fun. (e.g., 05 instead of 5)
      if row["Chromosome"].startswith("0"):
        row["Chromosome"] = row["Chromosome"][1:]
  

      coord = GenomeCoordinate(row["Chromosome"], None, row["ProteinStart"], row["ProteinStop"], row["ProteinID"])
      
      if coord in siteList:
        if coord not in dictOfDicts:
          genotypes = {}
          dictOfDicts[coord] = genotypes
        else:
          genotypes = dictOfDicts[coord]

        if haploid:
          genotypes[row["ID"]+"_A"] = row["DetectableAllele1"]
          genotypes[row["ID"]+"_B"] = row["DetectableAllele2"]
              
        else:
          genotypes[row["ID"]] = (row["DetectableAllele1"], row["DetectableAllele2"])
          

  fh.close()


  # make a synthetic diploid dataset from all pairs of haploid datasets
  if haploid:
    for coord in dictOfDicts:
      dictOfDicts[coord] = makeHaploidDiploid( dictOfDicts[coord] ) # take all pairs of haploid genotypes

  return dictOfDicts



def computeTheta(referenceGenotypes, coordinates, pops2samps, haploid2diploid=True):
  """
  This is a custom theta computation
  performed at the level of the locus.
  The locus is defined at the level of the chromosome;
  ie, considering all SAPs (translated into SNPs; or vice versa)
  but using the EXOME data.
  haplotypes that span the chromosome are at these sites are constructed
  and "theta" (Bi from Buckleton 2016; https://doi.org/10.1016/j.fsigen.2016.03.004)
  is computed and returned.
  This function is computed anew for each individual for every chromosome detected
  """
  
  if haploid2diploid==False:
    print("Not implemented yet..." , file=sys.stderr)
    sys.exit()

  popCounts = {} # population -> alleles -> counts (dictionary of dictionaries)
  popTots = {} # population -> diploidSamplesize (dictioary of string/int pairs)
  popHomozygosities = {} # this is M_i in Buckelton 2016
  # from this paper l is omitted as this computation is
  # done per locus on the fly.
  for pop, samples in pops2samps.items():

    # records, for a POPULATION, what the allele counts are
    hapCounts = defaultdict(int)

    nSamp = 0
    for sample in samples:
      haploidName1 = sample + "_A"
      haploidName2 = sample + "_B"

      hap1 = list()
      hap2 = list()
      for c in sorted(coordinates):
        dipName = haploidName1 + "\t" + haploidName2
        if dipName in referenceGenotypes[c]:
          (a1, a2) = referenceGenotypes[c][dipName]
        else:
          
          dipName = haploidName2 + "\t" + haploidName1
          if dipName not in referenceGenotypes: # some samples in our population LUT we do not have genotyping data for.
            continue
          (a1, a2) = referenceGenotypes[c][dipName]
          

        hap1.append(a1)
        hap2.append(a2)
        
      if len(hap1) == 0: # addresses samples present in the LUT and not found in the genotyping table
        continue

      nSamp += 2 # diploid...
      # convert the haplotypes to strings
      hap1 = ";".join(hap1)
      hap2 = ";".join(hap2)
      hapCounts[hap1] += 1
      hapCounts[hap2] += 1

      
    popCounts[pop] = hapCounts

    homozygosities = 0.0
    divisor = (nSamp)*( nSamp-1)
    popTots[pop] = nSamp
    # count is the number of times the haplotype was seen
    # the homozygosities is, given all pairs of individuals, the number of times
    # you see the same genotype. this is: count(count-1)/n(n-1)
    for count in hapCounts.values():
      homozygosities += (count*(count-1)) /divisor

    popHomozygosities[pop] = homozygosities # the haplotype homozygosity for population pop
    
  pops = list( pops2samps.keys() )
  betweenPopHomo=0. # This is Mb from Buckleton 2016
  # the between-population homozygosities
  
  for (pop1, pop2) in combinations(pops, 2):
    pop1Alleles = popCounts[pop1]
    pop2Alleles = popCounts[pop2]
    denom = float(popTots[pop1]*popTots[pop2])

    # this computes Mij (Buckleton)
    for a1 in pop1Alleles:
      # if not in pop2Alleles then c2 is 0, and numerator is 0 => no addition needed!
      if a1 in pop2Alleles:
        c1 = pop1Alleles[a1]
        c2 = pop2Alleles[a1]
        numerator = c1*c2
        
        betweenPopHomo += numerator / denom

  # in the paper they divide by |r|*(|r|-1)
  # which is the number of pairs of populations evaluated
  # BUT
  # they evaluate all pairs redundantly (as ordered pairs, CEU and CHB, then CHB and CEU)
  # our computation uses unordered pairs (CEU,CHB but not CHB,CEU)
  # we need to adjust by /2
  mbDivisor = (len(pops)*(len(pops)-1))/2
  Mb = betweenPopHomo / mbDivisor

  Bw = 0 # Bw from Buckleton
  for pop in pops:
    Mi = popHomozygosities[pop]
    if Mi > Mb: # take conditional sum (effectively truncate to remove negative between-population heterozygosities)
      Bw += (Mi-Mb) /(1-Mb)
    
    
  # Taken from Appendix A of Buckleton
  
  Bw /= len(pops)
    
  return max(Bw, ALT_MIN_THETA)

    
def getLogRmp(referenceGenotypes, originalCoordinatesToUse, diploidIndividuals, originalProteome, samps2pops, pops2samps, falsePositiveRate=0.01, numIters=1000):
  """
  This is the meat and potatoes of the RMP calculation. The interface is clunky
  I use a Monte Carlo estimate (numIter estimates of the RMP) and return the mean
  of these estimates. 
  Specifically, it returns the log RMP, the number of proteomic markers and the theta inferred

  The data provided to this function must be constrained to be within the same chromosome.
  This function starts by estimating a locus-specific theta (computeTheta)
  which in turn requires genotypes (referenceGenotypes), their coordinates (originalCoordinatesToUse)
  and a lookup table (dictionary) that associates a population to a list of individuals (pops2samps)
 
  Each iteration of the MC approach is as follows:
  The code iterates over all alleles within a chromosome (this is the locus)
  and drops alleles with rate falsePositiveRate
  and computes an exact RMP on this locus.
  In doing so it estimates the effect of drop-in (assuming independence of drop-in events)
  

  """

  
  if falsePositiveRate==0: #nothing to permute.
    numIters =1

  # get the number of proteomic markers
  nTot=0
  for c in originalCoordinatesToUse:
    p = originalProteome[c]
    
    if len(p)==1:
      nTot += 1
    else:
      nTot += 2

  
  theta = computeTheta(referenceGenotypes, originalCoordinatesToUse, pops2samps)
  # the RMPs based on imputation
  rmps = numpy.zeros(numIters)
  # and the RMPs based on the evidence
  rmpsEvidence = numpy.zeros(numIters)

  rmpsNaive = numpy.zeros(numIters)
  
  singleton = {}
  
  for i in range(numIters):
    coordinatesToUse = set() # the coordinates to evaluate in the MC simulation. subset of the original coordinates 
    proteome = dict()
    # a list, converted to string, used in the singleton design pattern
    # ensures that I only recompute the RMP when the protein diplotype is unique
    prots = []
    
    for c in originalCoordinatesToUse:
      p = originalProteome[c]
      
      # simulate drop-in
      # encoded as starting with an empty list and adding alleles
      # with probability 1- falsePositiveRate

      # this is made more complicated because the data are organized at the level of the segregating site
      # i.e., there can be 1 or 2 alleles detected per c
      
      if len(p)==1: # a single allele detected in proteome
        r = random.random()
        if r >= falsePositiveRate:
          coordinatesToUse.add(c)
          proteome[c]=p
          prots.append( str(c) + p[0].basecall )
      else: # individual is heterozygous in the proteome
        # p is of type list
        r1 = random.random()
        r2 = random.random()
        if r1 >= falsePositiveRate:
          if r2 >= falsePositiveRate: #keeping both proteomic alleles
            coordinatesToUse.add(c)
            prots.append( str(c) + p[0].basecall )
            prots.append( str(c) + p[1].basecall )
            proteome[c]=p
          else:# second allele in heterozygous site is dropped
            p2 = list()
            p2.append(p[0])
            prots.append( str(c) + p[0].basecall )
            coordinatesToUse.add(c)
            proteome[c]=p2
        elif r2 >= falsePositiveRate: # first allele...
          p2 = list()
          p2.append(p[1])
          prots.append( str(c) + p[1].basecall )
          coordinatesToUse.add(c)
          proteome[c]=p2

    # we've randomly dropped alleles from the proteome with probability falsePositiveRate
    protStr = ";".join(prots)
    if protStr not in singleton:
      # simple singleton design pattern
      r = getLogRmpExact(referenceGenotypes, list(coordinatesToUse), diploidIndividuals, proteome, falsePositiveRate=0, theta=theta)
      singleton[protStr] = r
    else:
      r = singleton[protStr]
      
    rmps[i] = exp(r[0])
    rmpsEvidence[i] = exp(r[2])
    rmpsNaive[i] = exp(r[3])
    
  meanRmpImputed = sum(rmps)/numIters # compute the expectation from the Monte Carlo simulation

  meanRmp = sum(rmpsEvidence)/numIters

  meanRmpNaive = sum(rmpsNaive)/numIters
  
  return(log(meanRmp), nTot, theta, log(meanRmpImputed), log(meanRmpNaive))
  

def getLogRmpExact(referenceGenotypes, coordinatesToUse, diploidIndividuals, proteome, falsePositiveRate=0.01, thetaCorrect=True, protBlacklist=None, theta=0.03):
  '''
  This is an exact computation 
  For now it is an ERROR to provide it a falsePositiveRate that is nonzero 
  In the pseudocode this is analogous to: rmpWithoutError
  save that the protein profile is now assumed to be 100% accurate.
  Several parameters are in the interface to this function that cannot and should
  not be used. 
  The false positive rate MUST be 0, and the thetaCorrection MUST be true.
  Future iterations may fix this, but for now the interface stands.
  
  In truth, error can be better incorporated. We know the number of mismatches
  for each exome sample (that is, the number and type of drop-in required 
  for each sample to be a match).
  We can use the binomial PMF to compute the RMP directly from this.
  HOWEVER
  there are nested sets that must be accounted for
  e.g., a perfect match in the exome space is also is consistent with 
  every possible subset of allelic dropout.
  handling this case is straight-forward
  but the other cases is sticky.
  '''  


  if falsePositiveRate!=0:
    print("No good. In the current implementation of *getLogRmpExact* an error rate of 0 is required. This may get relaxed later", file=sys.stderr)
    sys.exit()

  nTot=0 # total number of comparisons (haploids counting as 1, diploids as 2)
  maxMismatch=0
  mismatchDistibution = numpy.zeros( len(coordinatesToUse)*2 + 1, numpy.int64)

  
  for c in coordinatesToUse:
    p = proteome[c]
    if len(p) == 1:
      nTot += 1
    else:
      nTot += 2

  exactMatches = defaultdict(int) # the diplotypes consistent with the proteome (exact)
  overallAlleleCounts = defaultdict(int) # and the marginal allele counts (overall)
  haploidIDs = set()
  
  for dip in diploidIndividuals:
    numMismatches = 0# specifically, the number of false positives in the proteome needed to make these
    # proteome to be consistent with the diploid sample
    thisPair = dip.split("\t")
    haploidIDs.add(thisPair[0])
    haploidIDs.add(thisPair[1])
    
    for c in coordinatesToUse:
      if dip not in referenceGenotypes[c]:
        print("Should not be!", c, dip, sep="\n", file=sys.stderr)
        sys.exit()
        
      (a1,a2) = referenceGenotypes[c][dip]
      #p is either 1 or 2 proteome allele calls...
      p = proteome[c]
      if protBlacklist is not None: #
        if len(p) ==1:
          if p[0] in protBlacklist:
            continue
        else:
          if p[0] in protBlacklist:
            if p[1] in protBlacklist:
              continue
            newP = set()
            newP.add(p[1]) # being careful not to update p in proteome[c]
            p = newP
          elif p[1] in protBlacklist:
            newP = set()
            newP.add(p[0])
            p = newP
            
      if len(p) == 1: # haploid in proteome
        if p[0].basecall != a1 and p[0].basecall != a2: #and the allele is inconsistent
          numMismatches += 1
          if falsePositiveRate == 0:
            break
      else: # diploid in proteome
        if p[0].basecall != a1 and p[0].basecall != a2: # 0th basecall is inconsistent
          numMismatches += 1
          if falsePositiveRate == 0:
            break
        if p[1].basecall != a1 and p[1].basecall != a2: # 1st basecall is inconsistent
          numMismatches += 1
          if falsePositiveRate == 0: # 1 mismatch is sufficient to say it's not a match...
            break
          
    mismatchDistibution[ numMismatches ] += 1
    
    if numMismatches > maxMismatch:
      maxMismatch = numMismatches
      

    hap1  = list()
    hap2 =  list()
    for c in coordinatesToUse:
      (a1,a2) = referenceGenotypes[c][dip]

      hap1.append(a1)
      hap2.append(a2)


    hap1Str = ";".join(hap1)
    hap2Str = ";".join(hap2)

    # this individual was an exact matched
    if numMismatches==0:
    # ensure that the pair of haplotypes is consistently encoded
    # e.g., we don't know which allele is maternal and which allele is paternal
    # just that we have two haplotypes. That way, e.g., A/T and T/A get encoded as the same thing
    # (at the level of the haplotype; if we have AA (maternal) TT (paternal) we also want that to be equal to TT (maternal) and AA (paternal)
      if hap1Str< hap2Str:
        dipString = hap1Str + "/" + hap2Str
      else:
        dipString = hap2Str + "/" + hap1Str
        
      exactMatches[dipString] += 1
      
    overallAlleleCounts[ hap1Str ] += 1
    overallAlleleCounts[ hap2Str ] += 1

  nDip = float(len(ids)) # number of diploid samples...
  nOriginalHap = len(haploidIDs) # convert from haploid IDs into the original diploid sample size (without considering all pairs)

  if not thetaCorrect:

    print("This function is no longer supported...", file=sys.stderr)
    sys.exit(1)
    
    maxExact = max( exactMatches.values() ) # the most common diplotype
    mismatchDistibution[0] = maxExact # this is the correct procedure when the theta correction is neglected AND the error rate is 0.
    
    if protBlacklist != None or falsePositiveRate==0:
      if mismatchDistibution[0] < 5:
        return (log(5/(nDip+nDip)), nTot)
      return (log(mismatchDistibution[0]/(nDip+nDip)), nTot)

    print("Cannot be here", file=sys.stderr)
    sys.exit(1) # will never be here, but as a sanity check this is added...
  else:

    # the allele frequency (overallAlleleCounts[h1]/(2*n) is based on the nchoose2 combinations of
    # individuals. This impacts both the numerator and the denominator (equally). i.e.,
    # the allele frequency of a haplotype does not change when we consider all pairs of haplotypes
    # what *does* need to change is our estimate of a minimum allele frequency
    # we define this to be 5/2n, where n is the *original* diploid individual sample size
    altMin = 5.0/nOriginalHap
    # in case there are no exact matches, this is the RMP
    maxRmp = rmpHeterozygousWithTheta(theta, altMin, 1-altMin)
    
    sumRMPs = 0
    # evaluate the dictionary in sorted order (descending)
    for diplotype,count in exactMatches.items():#sorted(exactMatches.items(), reverse=True, key= lambda kv: kv[1]):
      (h1,h2) = diplotype.split("/") # twohaplotypes...
      if h1 == h2: # homozygous case
        if not h1 in overallAlleleCounts:
          print("Should never happen!", h1, diplotype, count, file=sys.stderr)
          sys.exit(1)
          
        pu = overallAlleleCounts[h1] / (nDip+nDip) # allele frequency
        rmp = rmpHomozygousWithTheta(theta, max(pu, altMin))

        
      else: # heterozygous case
        if not h1 in overallAlleleCounts or not h2 in overallAlleleCounts:
          print("Should never happen!", h1, h2, diplotype, count, file=sys.stderr)
          sys.exit(1)
          
        pu = overallAlleleCounts[h1] / (nDip+nDip) # allele frequency
        pv = overallAlleleCounts[h2] / (nDip+nDip) # allele frequency
        rmp = rmpHeterozygousWithTheta(theta, max(pu, altMin), max(pv, altMin))

      sumRMPs += rmp
      maxRmp = max(maxRmp, rmp)

    # apply the 5/2n rule when there are no exact matches
    if sumRMPs == 0:
      sumRMPs = maxRmp

    return(log(maxRmp), nTot, log(sumRMPs), log(    max(5, mismatchDistibution[0]) / nDip ) )
    
  rmp = None
  print("Cannot be here", file=sys.stderr)
  sys.exit(1)
  
  # this code is BROKEN
  # it doesn't account for the subsets mentioned in the docstring.
  # it also doesn't account for the distinct haplotypes that can be used to form the matching genotypes.
  for i in range( maxMismatch+1):
    # m is the number of inferred proteomic spectra that have i mismatches
    # i.e., the only way that the two samples are the same is if there was i mismatches
    # e.g., if i is 1, then we're talking about all the samples that are consistent with this
    # proteome sample if there was a single drop-in event.
    # m then is the number of times that this occurred in our genomic sample
    
    m = mismatchDistibution[i]
    if m == 0: # if fractionMatching==0 then the whole below becomes 0
      if rmp is None:
        rmp = 0
      continue

    # counts the number of ways to have 0 drop ins, then 1 drop in, then 2, ...
    numPossibilities = nChooseK(nTot, i)
    
    # true positive probability TIMES the number of true positives
    # and the false postive probability TIMES the number of false positives
    probMatch = ((1-falsePositiveRate)**(nTot-i)) * ((falsePositiveRate)**i)

    # and the number of "matching" profiles (with i drop ins)  converted into an empirical estimate of probability
    fractionMatchingProfiles = m/nDip
    
    # the RMP just considering profiles that matched with i inserted (false+) peptides
    thisRMP = numPossibilities * probMatch * fractionMatchingProfiles

    # take some care to ensure that None is returned when we no RMP computation has been performed.
    if rmp is None:
      rmp = thisRMP
    else: # the possibilties are mutually exclusive (exactly 1 drop in, exactly 2 drop ins, ... )
      rmp += thisRMP

  # apply the 5/2n rule
  if rmp < 5/nDip:
    return (log(5/nDip), nTot)
  
  return (log(rmp), nTot)



def parseSamps2Pops(filename, sampleCol, popCol):
  """
  This parses a TSV file of samples associated with their respective population and 
  returns 2 dictionaries; one associating a sample -> their pop (samps2pops)
  and another associating the pop -> the list of samples (pops2samps)
  e.g.,
  a sample ID (NA00015) is they to their corresponding population affiliation (dictionary)
  and it returns a population (TSI) to a list of IDs associated with that population
  
  The arguments needed is the filename to parse
  The file needs to be a TSV file needs with a header

  and columns are the *names* of the columns (args2 and 3)
  associated with the sample ID and the population IDs, respectively

  """

  samps2pops = {}
  pops2samps = defaultdict(list)
  with open(filename) as tsv:
    reader = csv.DictReader(tsv, delimiter="\t")
    for row in reader:
      if sampleCol not in row or popCol not in row:
        print(filename , " does not have the appropriate column names in it. expecting " , sampleCol , "\t" , popCol, file=sys.stderr)
        sys.exit()
      # associate a sample ID to their population
      # and a population to a list of sample IDs

      samps2pops[ row[sampleCol] ] = row[ popCol ]
      pops2samps[ row[ popCol ] ].append( row[ sampleCol])


  return( samps2pops, pops2samps)

if __name__ == "__main__":

  exampleUsage = "\nAn example usage of how to use this program is:\n" + sys.argv[0] + " -p peptideCalls.tsv -q sampsToPops.tsv -r 1KG.table.phased.tsv.gz > RMPs.tsv"

  parser = argparse.ArgumentParser(description="Let's compute some RMPs!\n" + exampleUsage)

  parser.add_argument('-r', '--reference_genotypes', dest='R', help="This is a *tsv* file (from vcf2table.py) of phased genotypes (nucleic acid space).", type=str, default='')
  parser.add_argument('-c', '--reference_consequences', dest='C', help="This is a *tsv* file (from gvp2null.py) of phased peptide consequences (protein space).", type=str, default='')
  parser.add_argument('-p', '--peptide_hits', dest='P', help="This is a file of peptides detected augmented with their expected alleles in the nucleic acid space", type=str, default='')
  parser.add_argument('-q', '--population_to_samples', dest='Q', help="This is a file gives the population IDs of each sample specified in the file R", type=str, default="sampsToPops.tsv")
  parser.add_argument('-s', '--sample', dest='S', help="This restrains the analysis to just the sample id S", type=str, default='')
  parser.add_argument('-S', '--seed', dest='D', help="This sets the seed on the random number generator; helps with reproducibility", type=int, default=1)
  parser.add_argument('-l', '--peptide_lut', dest='L', help="A lookup-table to give to relate alleles to loci in peptide coordinates", type=str, default="peptideLUT.tsv")
  parser.add_argument('-e', '--error_rate', dest='E', help="Error rate(s) (i.e., drop-in rates) used in the RMP calculation.",type=float, default=[0, 0.01, 0.02, 0.03, 0.04, 0.05], nargs="+")
  
  results = parser.parse_known_args(sys.argv[1:])[0]
  args = parser.parse_known_args(sys.argv[1:])[1]

  genotypeTableFile = results.R
  peptideTableFile = results.C

  errs = results.E

  if min(errs) < 0:
    print("Problem with error rates... at least one is negative", errs, sep="\n", file=sys.stderr)
    sys.exit(1)
  if (max(errs) > 1):
    print("Problem with error rates... at least one is more than 1", errs, sep="\n", file=sys.stderr)
    sys.exit(1)
  
  if os.path.isfile(genotypeTableFile) and os.path.isfile(peptideTableFile):
    print("I need either a genotype table file (-r) or a protein consequences file (-c); not both!", file=sys.stderr)
    parser.print_help()
    sys.exit(1)
  
  if not os.path.isfile(genotypeTableFile) and not os.path.isfile(peptideTableFile):
    print("I need either a genotype table file (-r) or a protein consequences file (-c); you gave me neither", file=sys.stderr)
    parser.print_help()
    sys.exit(1)

  mode = NUCLEIC_ACID
  if os.path.isfile(peptideTableFile):
    mode = PROTEIN
    if not os.path.isfile(results.L):
      print("I need a peptide lookup table to work in peptide mode!", "Cannot find: ", results.L, file=sys.stderr, sep="\n")
      parser.print_help()
      sys.exit(1)

  proteomeHitsFile = results.P

  if not os.path.isfile(proteomeHitsFile):
    print("There is no file: " , proteomeHitsFile , file=sys.stderr, sep="\n")
    parser.print_help()
    sys.exit(1)

  prots = organizeProt( proteomeHitsFile, mode, results.L)
  
  samples = sorted(prots.keys() )
  
  sampsToPopsFile = results.Q
  if not os.path.isfile(sampsToPopsFile):
    print("There is no file: " , sampsToPopsFile , file=sys.stderr, sep="\n")
    parser.print_help()
    sys.exit(1)

  (samps2pops, pops2samps) = parseSamps2Pops(sampsToPopsFile, "Individual ID", "Population")
  
  if results.S != "":
    if results.S not in prots:
      print("There is no sample id:", results.S, file=sys.stderr, sep="\n", end="\n\n")
      print("The sample ids are:", "\n".join(sorted(prots.keys()) ), file=sys.stderr, sep="\n", end="\n\n")
      parser.print_help()
      sys.exit(1)
               
    samples = [results.S]
    
  
  print("SampleID", "Chromosome", "rmp", "nProtUsed", "theta", "imputed_rmp", "naive_rmp", "ErrorRate", sep="\t")
  
  numpy.random.seed(results.D) 
  
  for sampleID in samples: # look at each sample's proteome hits
    d = prots[sampleID]

    if mode == NUCLEIC_ACID:
      gts = loadGenotypes(genotypeTableFile, d, samps2pops, True) # and load genotypes (Haploid==True) from the 1000 genomes project
    else:
      gts = loadGenotypes(peptideTableFile, d, samps2pops, True, PROTEIN) # and load genotypes (Haploid==True) from the 1000 genomes project

    chromosomes = defaultdict(set)
    ids = set()
    smallestMarginalRmps= {}
    lowerboundRmp = {}
    for coord in gts:
      chromosomes[ coord.chromosome ].add(coord)
      for i in gts[coord]:
        ids.add(i)

    numDip = len(ids)
    
    chroms = sorted(chromosomes.keys())
    for chrom, markers in chromosomes.items():
      for err in errs:
        logRMP = getLogRmp(gts, markers, ids, d, samps2pops, pops2samps, err)
      
        print(sampleID, chrom, exp(logRMP[0]), logRMP[1], logRMP[2], exp(logRMP[3]), exp(logRMP[4]), err, sep="\t")

      
    
    

