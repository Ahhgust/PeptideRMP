#!/usr/bin/python3
# Written by August Woerner
# This program is designed to take in
# a peptide alleles, and generate
# a table of which individuals have these alleles.
# It was designed to work on Ensembl annotations
# in conjunction with the output of bcftools csq
#
# This table serves as a substrate for down-stream
# popgen statistics (RMPs, Fst).
# To perform the calculation I need to know
# the peptide allele(s): (SFSLFLSDGQR)
# which protein it came from: (ENSP00000257198)
# the reference coordinates 1-based, inclusive: (91      101) 
# and the chromosome: (I may ditch this later, as it's derivable from the protein ID)
#
# This program will group the alleles into loci, and infer the predicted peptide
# sequence for each haploid individual (as presented by bcftools csq)
#
# From this I will generate a table of peptide calls
# each row will represent a phased diploid individual's predicted peptide sequence
# Thus the two columns of peptide sequencing calls will, say, represent the maternal (first column)
# and paternal (subsequent column) peptide calls under the assumption that all of the genes
# are turned ON in this individual.
#
# CUrrent limitations:
# ONly autosomal calls are supported
# Only BCFtools csq is supported
# Only Ensembl annotations are supported



import sys
import gzip
import csv
import os.path
import glob
from collections import namedtuple, defaultdict
import argparse



#protFile = "Homo_sapiens.GRCh37.pep.all.fa.gz"
#gvpPanelFile = "gvpPanel.tsv"
#csqFile = 'ProteinPredictions/AllSimpleConsequences.txt.gz'

protCol = "ensembl_protein_id"
peptideCol = 'peptide_seq'
peptideStart = 'start'
peptideStop = 'end'
genomeCoordColumn = 'hg38start'



# lazy person's objects...
Allele = namedtuple('Allele', ['ID', 'protStart', 'protStop', 'seq', 'originalSeq'])
Variant = namedtuple("Variant", ['fromAllele', 'fromPos', 'toAllele', 'toPos'])


CigarOp = namedtuple("CigarOp", ["op", "size"])

def getPair(v):
  """
  Parses a variant (9338V)
  and splits it into its numeric and alphabetic bits...
  """

  i=0
  vlen = len(v)
  num = 0
  
  while i < vlen:
    if v[i] >= '0' and v[i] <= '9':
      num = (num*10) + int(v[i])
      i += 1
    else:
      break

  # some variants have a trailing * to indicate the end of protein sequence...
  #if v[-1] == '*' and i < vlen - 1:
    #return(num, v[i:-1])
  return(num, v[i:])


def parseVariant(s):
  '''
  This parses a variation column from bcftools csq
  eg, 9338V>9338T or
  (meaning at 1-based index position 9338, the V gets changed to a T)
  or
  9338VTT>9338V
  (meaning TT deletion)

  a Variant namedtuple is returned.
  the coordinates in the Variant (fromPos, toPos) are in zero-based
  half open coordinates (fromPos is 0-based, toPos is 1-based).
  in general protein[fromPos:toPos] should refer to the correct substring.
  '''
  
  (f, t) = s.split(">")
  (fromPos, fromVar) = getPair(f)
  (toPos, toVar) = getPair(t)

  # Wrong; the toVar position is the substring position in the individual, not the reference
  # e.g., it shifts after indels
  # we want the position of the variant in the reference in Bed coordinates
  #  return Variant(fromVar, fromPos-1, toVar, toPos)
  diffLen = len(fromVar) - len(toVar)
  if diffLen > 0:
    return Variant(fromVar, fromPos-1, toVar, fromPos + diffLen)

  return Variant(fromVar, fromPos-1, toVar, fromPos)
  

def getHaplotypes(seqsDict, trans2prot, consequencesFileOrGlob):
  """
  This parses the results of bcftools csq 
  (consequences file)
  and it takes a dictionary of ids -> protein sequences (reference) (seqsDict)
  and a LUT that translates transcript IDs to protein ids (dictionary)
  and it returns a dictionary that associated a haploid ID
  with the variants associated with the haploid sequence for all
  proteins in seqsDict
  """

  # may be 1 file, or a glob 'file*txt' for many files...
  csqFiles = [consequencesFileOrGlob]

  if not os.path.isfile(consequencesFileOrGlob):
    csqFiles = glob.glob(consequencesFileOrGlob)
    if not csqFiles:
      print("The consequences file(s):",  consequencesFileOrGlob, "isn't a file and it also isn't a unix-style pattern for files... what to do!", file=sys.stderr, sep="\n")
      sys.exit(1)
      

  d = {}
  ids = set()

  for consequencesFile in csqFiles:
  
    if consequencesFile.endswith('.gz'):
      fh = gzip.open(consequencesFile, "rt")
    else:
      fh = open(consequencesFile, 'r')

  
    for line in fh:
      if line.startswith("#") or line.startswith("LOG"):
        continue
      s = line.rstrip().split("\t")
      geneInfo = s[-1].split("|")
      if len(geneInfo) < 3:
        print(line, geneInfo)
        sys.exit()
      transcript = geneInfo[2]
      # this is one of the transcripts we're after
      ids.add(s[1])
      
      if transcript in trans2prot:
      
        prot = trans2prot[transcript]
        who = s[1] + "_" + s[2] # sampleID _ chromosome (maternal or paternal)
        what = parseVariant(geneInfo[5])

        if geneInfo[0] == 'frameshift' and geneInfo[5][-1] != '*':
          geneInfo[5] += "*" # concatenate a separate *, indicating that the protein is "done"
        # on occassion the * is missing from the annotation with frameshift mutations, apparently
        # the rules are csq looks for a stop codon; if it finds one that's not much farther than
        # the original stop codon's position, it adds a '*'. Otherwise the sequence just stops.
        # this adds the "*" annotation back in. 
      
        # bcftools csq encodes an event as, say, *missense if it would ordinarly be encoded as a
        # missense variant, BUT there's an upstream variation (e.g., frameshift) that causes
        # this to be a null allele.
        # commented b/c this should already be apparent-- if it's a missense b/c of a frameshift, just
        # use the frameshift annotation.


          
        if what.fromPos < 0: # (effectively) synonymous variation... (has no index in the protein)
          continue

        if prot not in d:
          d[prot] = defaultdict(list)

        
        innerDict = d[prot]
        # sometimes the same variant has multiple characterizations
        # when this happens its the same protein coordinates that occur twice in a row. 
        # let's ensure that's not happening either (just add the variant once)
        if len(innerDict[who])== 0 or (innerDict[who][-1].toPos != what.toPos and innerDict[who][-1].fromPos != what.fromPos):
          innerDict[who].append(what)

    fh.close()


  
  for locus in d:
    innerDict = d[locus]
    # some protein coordinates are on the negative strand, so
    # from the perspective of the genome the difference encoding
    # is occurring on reverse sorted order
    # let's make sure that everything is sorted (ascending)
    for whatList in innerDict.values():
      whatList.sort(key=lambda variant: (variant.fromPos, variant.toPos) )
      
  
  # a list of sample IDs (deterministically ordered)
  idsOrdered = list(ids)
  idsOrdered.sort() # we need a definitive list of who was genotyped. As a proxy, we take this to be the
  # union of all individuals in the difference encodings.
  # strictly this may miss any individual that is homozygous reference for everything... but that's
  # quite unlikely.
  
  return d, idsOrdered

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

def validateCigar(cigarOps, s):
  slen = len(s)
  tot=0
  for c in cigarOps:
    if c.op != 'D':
      tot += c.size
      
  if tot != slen:
    print("Should never happen\n", tot, slen, cigarOps, s, sep="\n", file=sys.stderr)
    return False
  return True

def cigarops2string(cigarOps):
  """
  Takes in a list of cigar tuples (c.op in the set MID, c.size is an int)
  and returns a cigar string (all concatenated together)
  """
  # common case
  if len(cigarOps)== 1:
    return str(cigarOps[0].size) + cigarOps[0].op
  # general
  return "".join( [str(c.size) + c.op for c in cigarOps] )

def diffs2seq(diffs, seq, cigar=False):
  """
  This takes a sequence of difference encodings (variants)
  and a reference sequence
  and returns a 2-pl
  the first is a string (with isoleucines converted to leucines)
  the second is a list (without the conversion)
  the list maintains the genomic indexing (really proteomic indexing)
  so that slices may be taken at a given index and the correct sequence is returned.

  This is complicated by indels, which have a genomic index but that index
  may not be unique. the slice approach is approximate in that respect
  while the string returned is just that. a string.
  it's the string that the different alleles are searched in,
  and the slice notation that's used when we fail to find a matching allele
  """

  newSeq = list(seq)
  cigarOps = []
  truncateTo = -1

  # keep track of deletions and insertions
  numDeleted=0
  numInserted=0

  # edge case; no diffs, and there's a stop-codon
  # the other case to worry about is that there are diffs
  # in which case the stop-codon can be altered (or not)
  # happens with the polymorphic pseudogene class
  if not diffs:
    for i in range(len(newSeq)):
      if newSeq[i] == '*':
        truncateTo = i
        break
  
  for d in diffs:

    # coordinates wrt to the reference sequence
    x = d.fromPos
    # check all positions inthe protein PRIOR to this point
    # and see if there's a stop-codon
    for i in range(x):
      if newSeq[i] == '*':
        truncateTo = i
        break
      
    if truncateTo >= 0:
      break
    # STOP codon...
    if d.toAllele is None:
     # return (isoleucine2leucine( "".join(newSeq[:x])), newSeq[:x])
      truncateTo = x
      break
    
    # proteomic index
    y = x + len(d.fromAllele)
    # substring position within the from allele
    j = 0
    
    # overstrike the bases in the range
    for i in range(x,y):
      # deletion wrt reference
      if j < len(d.toAllele):
        if d.toAllele[j] == '*': # end of protein
          #return (isoleucine2leucine( "".join(newSeq[:i])), newSeq[:i])
          truncateTo = i
          break
          
        if i < len(newSeq):
          newSeq[i] = d.toAllele[j]
        else:
          y = i
          break # cause the if statement to treat this as a concatenation...
      else: # case of deletion...
        if i < len(newSeq):
          newSeq[i] = ''
          numDeleted += 1
          if i > 0 and newSeq[i-1] == '': # update the size of the deletion
            cigarOps[-1] = CigarOp('D', cigarOps[-1].size + 1)
          else:
            # +1 for the pre-increment of numDeleted
            nMatches = getPreviousMatchlen(cigarOps, i - numDeleted + numInserted + 1)
            if nMatches > 0:
              cigarOps.append( CigarOp('M', nMatches) )
              
            cigarOps.append( CigarOp('D', 1) )
          
      j += 1

    if truncateTo != -1:
      break
    # insertion event.
    # if the toAllele is longer, concatenate it (this is an insertion wrt the reference)
    if j < len(d.toAllele):

      # +1 b/c we need to include that ith base as matching...
      # but +1 is implicit b/c we're on index y, but we need to refer to index y-1 
      nMatches = getPreviousMatchlen(cigarOps, y - numDeleted + numInserted)
      if nMatches > 0:
        cigarOps.append( CigarOp('M', nMatches) )
      
      if d.toAllele[-1] == '*': # don't need the truncation marked
        newAllele = d.toAllele[j:-1]
        cigarOps.append( CigarOp("I", len(newAllele)))
        newSeq[y-1] += newAllele
        truncateTo = y
        break  
              
      else:
        newAllele = d.toAllele[j:]
        cigarOps.append( CigarOp("I", len(newAllele)))
        newSeq[y-1] += newAllele
        numInserted = cigarOps[-1].size

  terminalDel = 0
  if truncateTo >= 0:
    terminalDel = len(newSeq) - truncateTo # only works b/c no nucletides have been inserted/deleted in the sequence AFTER the truncation
    newSeq = newSeq[ : truncateTo ] 

  i2l = isoleucine2leucine("".join(newSeq) )
  nMatches = getPreviousMatchlen(cigarOps, len(i2l))

  if nMatches > 0:
    cigarOps.append( CigarOp('M', nMatches) )

  if terminalDel:
    cigarOps.append( CigarOp('D', terminalDel) )

  if not validateCigar(cigarOps, i2l):
    print(diffs, truncateTo, newAllele, sep="\n", file=sys.stderr)
    exit(1)

    
  if cigar:
    return (i2l, newSeq, cigarOps)
  return (i2l, newSeq)




def seqToProteome(s):
  """
  Takes in a list or string 
  and it returns a string
  that has been (re) digested by typsin
  and the Isoleucines have been converted into Leucines
  """
  
  if type(s) == str:
    outList = s.split("")
  elif type(s) == list: # convoluted. s is a list of strings. I want a list of chars...
    outList = "".join(s).split("")
  else:
    print("Unsupported object " , s , file=sys.stderr)
    return ""
  
  sLen = len(outList)
  for i in range(sLen):
    c = outList[i]
    if willDigest(c):
      return "".join(outList[:i])
    elif c == 'I': # leucine and isoleucine are mass-identical. convert!
      outList[i] = 'L'
      
  return "".join(outList)



def makeLoci(alleles):
  '''
  This takes the alleles associated with a given protein
  and them returns a list of Loci
  
  The allele locations (loci) are evaluated and partitioned
  s.t. each locus is associated with at least 2 alleles
  This includes 1+ alleles that are sought
  as well as an additional placeholder (Null) allele
  that represents all alleles that were not sought.
  '''

  nAll = len(alleles)
  j = i = 0
  
  listOfTups = []
  while (i < nAll):
    j = i + 1
    while j < nAll:
      # different alleles of the same locus
      # ie loci overlap.
      if alleles[i].protStart < alleles[j].protStop and \
         alleles[i].protStop > alleles[j].protStart:
        j += 1
      # one allele is strictly before the other
      elif alleles[i].protStop <= alleles[j].protStart:
        break
      else:
        print(alleles[i] , "is not distinct from" , alleles[j] , "The allelic regions need to be distinct!", file=sys.stderr, sep="\n")
        sys.exit()

    x = alleles[i].protStart
    y = alleles[i].protStop
    for z in alleles[(i+1):j]:
      if z.protStart < x:
        x = z.protStart
      if z.protStop > y:
        y = z.protStop
      
    newList = [ Allele( alleles[i][0], x, y, "null", "null")  ]
    # the reference allele can be present multiple times
    # making this into a set helps out with that.
    newList.extend(set(alleles[i:j]))
    #print(newList[0])
    # do a sanity check.
    # check all pairs ensuring that either the first coordinate matches or the second...
    nLoci = len(newList)
    
    for x in range(1, nLoci):
      locus = newList[x]
      for y in range(x+1, nLoci):
        if locus.protStart != newList[y].protStart and \
           locus.protStop != newList[y].protStop:
          print("There may be a clustering problem-- generally alleles in a locus share either a start or stop position", \
                "This is not the case of",\
                locus , "vs" , \
                newList[y],\
                file=sys.stderr, sep="\n")
    
    listOfTups.append( tuple( newList ) )
    i = j

  return listOfTups



def getDictOfSeqs(protFilename, protDict):
  '''
  This parses a fasta file of peptide sequences as downloaded
  from Ensembl. 
  e.g., file:
  Homo_sapiens.GRCh37.pep.all.fa.gz
  This is protFilename
  As an optimization, only protein IDs in the protDictionary are considered
  The protDictionary may be None (in which case all proteins are returned)

  Embedded in the fasta file header is the name of the protein (ENSP)
  transcript (ENST) and gene (ENSG)
  This is courtesy of Ensembl.

  I don't care about the gene name, but translating between ENST and ENSP (2nd return value)
  and vice versa (3rd return) while the first return value is a dict that: ENSP->protein sequence
  '''  
  d = {}

  # trans ID to PROT id
  translator = {}
  # vice versa
  protein2trans = {}
  
  if protFilename.endswith(".gz"):
    f = gzip.open(protFilename, "rt")
  else:
    f = open(protFilename, 'r')

  prot = ''
  trans = ''
  seq = ''
  for s in f:
    # boiler plate parse FASTA + recover the gene/protein/transcript names...
    if s.startswith('>'):

      if seq != "" and trans != '':
        d[prot] = seq
        translator[trans] = prot
        protein2trans[prot] = trans
        
      sep = s[1:].split(' ')
      prot = sep[0].split(".")[0]
      seq = trans = ''
      if protDict is not None and prot not in protDict:
        continue # trans == '' if this is not a protein we want

      for t in sep:
        if t.startswith("transcript:"):
          trans = t[11:].split(".")[0]
          break

      if trans == "":
        print("Shouldn't happen: " , sep, file=sys.stderr, sep="\n")

    else:
      seq += s.rstrip()
      

  f.close()

  # grab tthe trailing value... if it is one we want
  if seq != "" and trans != '':
    d[prot] = seq
    translator[trans] = prot
    protein2trans[prot] = trans
  # make sure every gene we asked for, we got!
  errs=0
  if protDict is not None:
    for prot in protDict:
      if prot not in d:
        errs+=1
        print("Missing " , prot , sep="\n", file=sys.stderr)

    if errs:
      sys.exit(1)
  
  return d, translator, protein2trans
  


def willDigest(char):
  """
  Overload this to try out other digests.
  Right now it just takes single character as an argument...
  and tells you if it's digest-able.
  Currently only a trypsin digest is encoded.
  """
  
  if char == 'R' or char == 'K':
    return True
  return False


def digest(s, asList=False):
  """
  Takes in a string or a list
  and returns the protein trypsin digest of s
  it also converts Isoleucines to Leucines
  
  If asList is True then the first peptide from the digest is returned
  otherwise a list of all peptide products is returned.
  """
  
  if type(s) == list:
    s = "".join(s)

  listOfLists = []
  l = []

  for c in s:
    if c == 'I':
      l.append("L")
    else:
      l.append(c)

    if willDigest(c):
      if asList == False:
        break
      else:
        listOfLists.append( "".join(l) )
        l = []

  if len(l):
    listOfLists.append(l)

  if asList == False:
    return "".join(l)

  return listOfLists



def getRefAllele(allele, protein, basePrior=False):
  '''
  This routine takes in an allele object and a protein string
  and returns the reference sequence associated with that allele.
  the reference sequence is returned as a list

  If basePrior is True then the 5' (direction of N terminus)
  amino acid is also returned 
  This is useful to return as this character should be digestable...

  '''
  
  if allele.protStop > len(protein):
    print(allele, " exceeds the protein length of " , len(protein), \
          "that doesn't sound right!", sep="\n", file=sys.stderr)
    sys.exit(1)

    
  if allele.protStop != len(protein):
    lastChar = protein[ allele.protStop - 1]
    if not willDigest(lastChar):
      print("The reference allele of " ,\
            allele, protein[allele.protStart:allele.protStop], \
            "does not end in an R or K site. This may or may not be a problem...", \
            sep="\n", file=sys.stderr)

  if allele.protStart > 0:
    firstChar = protein[ allele.protStart - 1]
    if not willDigest(firstChar):
      print("The reference allele of " ,\
            allele, protein[allele.protStart:allele.protStop], \
            "does not have an R or K site prior to its start. This may or may not be a problem...", \
            sep="\n", file=sys.stderr)
    

  if basePrior and allele.protStart > 0:
    return list(protein[(allele.protStart-1):allele.protStop])

  # doubtful if this can occur; first AA is encoded as Methionine in Ensembl
  # this may not be actually part of the peptide.
  return list(protein[allele.protStart:allele.protStop])


def getConsistentAllele(alleles, a):
  '''
  This takes a list of alleles where alleles[0] is the null allele
  alleles are defined at a particular peptide spanning 1+ SAPs
  and it takes the full haploid protein sequence (a)
  and it checks the number of allels that match a 
  This is done by digesting the protein
  and then looking for an exact match 
  in the set of alleles 
  
  If the peptide sequences sought are correctly identified then
  there should be exactly 0 (NO allele detected)
  or exactly 1 (1 consistent allele detected) hit
  If multiple hits are found this implies that the peptide sequence sought
  is multi-copy (within the protein), which is not currently permitted
  '''
  
  matches = []

  substrings =  digest(a, asList=True) 
  
  for putativeMatch in alleles[1:]:
# this was a coordinate-based approach. I don't think it's the right solution...
    # the substring positions
#    sStart = 0
 #   sStop = stop - start
    
  #  if start < putativeMatch.protStart:
   #   sStart = putativeMatch.protStart - start
    #if stop > putativeMatch.protStart:
     # sStop = putativeMatch.protStart - start

#    substrings = digest(a[sStart:sStop], asList=True)

    if putativeMatch.seq in substrings:
      matches.append( putativeMatch )



  if len(matches) == 0:
    return alleles[0] # null allele
  elif len(matches) == 1:
    return matches[0]
  
  print("Not good. The peptide: ", "".join(a), "matches these alleles: " , matches ,\
        "As in more than one. Not good!", file=sys.stderr, sep="\n")
  return None
  


def getPeptideHaps(proteinHaps, alleles, refSeq, proteinSeq):
  '''
  This takes in a dictionary of protein haplotypes (key => sampleID, value => list of sequence differences from the reference sequence for sampleID
  a locus, which is a tuple of alleles where the 0th tuple is for the null allele which spans the sequence coordinates of all alleles
  and the peptide sequence of the null allele as a list
  '''
  
  d = {}
  nullAllele = alleles[0] # this locus *spans* all alleles at this locus
  refAsString = "".join(refSeq)
  
  uniqueHaps = {}

  errs = 0
  for sample, diffs in proteinHaps.items():

    locus = []
    n = len(diffs)

    (haplotype, hapAsList) = diffs2seq(diffs, proteinSeq)

    matches = []
    for a in alleles[1:]:

      # both sequences have had Isoleucines converted to Leucines
      f = haplotype.find(a.seq)
      if f != -1:
        if f == 0 or willDigest(haplotype[f-1]):
          matches.append(a)

    if len(matches) == 0:
      allele = alleles[0]
    elif len(matches) == 1:
      allele = matches[0]
    else:
      errs += 1
      if errs == 1:
        print("Should never happen.", sample, "The peptides at locus", alleles[1:] , "are problematic-- alleles within a locus should be mutually exclusive, these are not!E.g.,", haplotype, sep="\n", file=sys.stderr)
      # called as a null allele...
      # but added a separate tag
      badAllele= Allele( alleles[0][0], alleles[0][1], alleles[0][2], "notdistinct", alleles[0][4])
      allele = badAllele

    x = alleles[0].protStart
    y = alleles[0].protStop

    # needed when grabbing the first peptide!
    if x < 1:
      x = 1
    
    variant = ""
    if x >= len(hapAsList): # variant is strictly AFTER any sequence in this peptide
      variant = ""
    elif y >= len(hapAsList): # variant begins before the peptide ends, but ends AFTER it...
      variant = "".join(hapAsList[(x-1):])
    else:
      variant = "".join(hapAsList[(x-1):y]) # common case; variant is a substring 

    d[sample] = (variant, allele)
    continue
    # below is code used to try and get the same information using positional information on
    # the variation. In the end this approach is flawed because the SNP/indel positions are not
    # unique and may be arbitrary. I thought I could approximate a solution based on the SNPs/
    # indels within the peptide region... but I was wrong.
    diffJustPrior = None
    
    for a in range(n):
      diff = diffs[a]

      if diff.toAllele is None:
        locus = None
        break
      
      i = diff.fromPos - nullAllele.protStart
      j = nullAllele.protStop - diff.fromPos
      
      f = diff.fromAllele
      t = diff.toAllele

      m = max(len(f), len(t))
      
      if i >= 0 and j >= 0: # variation is strictly contained within the peptide
        if len(locus) == 0:
          locus = refSeq[:] # deep copy, yo

        j = i + m
        y= 0
        for x in range(i,j):
          if x == len(locus):
            locus.append("")
            
          if i == j-1:
            locus[x] = t[y:]
          else:
            locus[x] = t[y:y] # slice notation gives a "" if I go beyond the end of the string...
          y += 1
      elif j < 0:
        break
      # overlaps 5' of peptide
      elif diff.toPos >= nullAllele.protStart and diff.fromPos <= nullAllele.protStart:
        # handle the case of insertion and missense:
        if diff.fromPos == diff.toPos -1:
          if not willDigest(diff.toAllele[-1]):
            print("whoops; lost it")
            
        print(diff)
        print(nullAllele)
        sys.exit()
      
      diffJustPrior = diff

    # upstream variant obliterated this protein sequence
    if locus is None:
      d[sample] = ("stop", alleles[0])
    elif len(locus)==0:
      if refAsString not in uniqueHaps:
        uniqueHaps[ refAsString ] = getConsistentAllele(alleles, refSeq)
        
      d[sample] = (refAsString, uniqueHaps[ refAsString ] )
    else:
      lString = "".join(locus)
      if lString not in uniqueHaps:
        uniqueHaps[ lString ] = getConsistentAllele(alleles, locus)
      d[sample] = ("".join(locus), uniqueHaps[ lString ] )      

    
  return d


def uniquify(s, delimiter=";"):
  """
  Helper function
  It simplifies a delimited string with redundant fields
  if one unique field exist, a string is returned
  otherwise a list of the unique elements from the string is returned
  """

  s = s.split(delimiter)
  s = list(set(s))
  if len(s)==1:
    return s[0]
  
  return s


def getPreviousMatchlen(cigarOps, currOffset):
  """
  This takes in a vector of cigar operations (reduced instruction set; MID only)
  and returns the length of the match needed to consume bases 0..currOffset in the string
  """

  if not cigarOps:
    return currOffset

  off = 0
  # recompute the current substring position in the read 
  for cig in cigarOps:
    if cig.op != 'D':
      off += cig.size

  # and return the amount needed to be consumed...
  return currOffset - off

  

def dumpProts(protFile, consequencesFile):

  
  (seqsDict, trans2prot, prot2trans) = getDictOfSeqs(protFile, None)
  (haps, ids) = getHaplotypes(seqsDict, trans2prot, consequencesFile)

  fh = gzip.open("dumpLUT.tsv.gz", "wt")
  fh.write("ProteinID\tSampleID\tCigar\n")

  
    
  for protID, refSeq in seqsDict.items():

    refSeq = isoleucine2leucine(refSeq)
    
    # associate a protein haplotype with who has it
    proteinHaplotypes = defaultdict(list)
    # and a protein haplotype with the difference encoding
    proteinDiffs = dict()

    if protID in haps:
      protHaps = haps[ protID ]
      for who in ids:
        who_1 = who + "_1"
        who_2 = who + "_2"
        
        if who_1 in protHaps:
          diffs = protHaps[who_1]
          (seq, listy, ciggy) = diffs2seq( diffs , refSeq, cigar=True)
        else:
          #seq = refSeq
          #ciggy = [ CigarOp("M", len(seq)) ]
          (seq, listy, ciggy) = diffs2seq( [] , refSeq, cigar=True)
          
        # re-add the stop-codons back in.
        # these are present iff the original Ensembl sequence had them
        # e.g., in the case of polymorphic pseudogenes
        offset = seq.find("*")
        if offset < 0:
          proteinHaplotypes[seq].append( (who_1, ciggy) )
        else:
          ciggy = [ CigarOp("M", offset) ]
          proteinHaplotypes[seq[ : offset ] ].append( (who_1, ciggy) )


        
        # and repeat for the second allele...
        if who_2 in protHaps:
          diffs = protHaps[who_2]
          (seq, listy, ciggy) = diffs2seq( diffs , refSeq, cigar=True)
          
        else:
          #seq = refSeq
          #ciggy = [ CigarOp("M", len(seq)) ]
          (seq, listy, ciggy) = diffs2seq( [] , refSeq, cigar=True)
          
        offset = seq.find("*")
        if offset < 0:
          proteinHaplotypes[seq].append( (who_2, ciggy) )
        else:
          ciggy = [ CigarOp("M", len(seq)) ]
          proteinHaplotypes[seq[ : offset ] ].append( (who_2, ciggy) )

          

    else: # need to just print out the ref sequence...
      bigList = [] # which is associated with everyone...
      
      (seq, listy, ciggy) = diffs2seq( [] , refSeq, cigar=True)
      for who in ids:
        bigList.append( (who + "_1", ciggy) )
        bigList.append( (who + "_2", ciggy) )

      proteinHaplotypes[seq] = bigList


    i = 0
    for seq, whos in proteinHaplotypes.items():
      if seq == "":
        continue
      i += 1
      print(">", protID , "_" , i, sep="")
      print(seq)
      
      for who in whos:
        fh.write(protID + "_" + str(i) + "\t" + who[0] + "\t" + cigarops2string(who[1]) +  "\n")

  fh.close()

  

def main(argv):
  parser = argparse.ArgumentParser(description="Let's get some peptides!")

  parser.add_argument('-g', '--gvp_panel', dest='G', help="The GVP panel file", type=str, default="gvpPanel.tsv")
  parser.add_argument('-c', '--csq_file', dest='C', help="The output from bcftools csq (gzipped and grepped to just included protein-changing variation)", type=str, default="AllSimpleConsequences.txt.gz")
  parser.add_argument('-p', '--protein_references', dest='P', help="A \"pep\" file from Ensembl. e.g., \"Homo_sapiens.GRCh37.pep.all.fa.gz\". Must correspond to -c (as in the Ensembl protein coordinates must correspond to the consequences file!) ", type=str, default="Homo_sapiens.GRCh37.pep.all.fa.gz")
  parser.add_argument('-i', '--ignore_gvps', dest='I', help="A comma separated list of GVP *alleles* that are to be ignored/dropped/omitted from the output", type=str, default="")
  parser.add_argument('-l', '--lut', dest='L', help="The lookup table to convert the coordinates of the allele(s) to that of the locus", type=str, default="peptideLUT.tsv")
  parser.add_argument('-d', '--dump', dest='D', help="Ignores the gvp_panel and just dumps the fasta sequence of each haploid individual", action='store_true', default=False)

  
  results = parser.parse_known_args(argv[1:])[0]
  args = parser.parse_known_args(argv[1:])[1]
  if args:
    print("Extra arguments detected", args , sep="\n", file=sys.stderr)
    return 1
  
  gvpPanelFile = results.G
  csqFile = results.C
  protFile = results.P
  
  ignoredAlleles = {}
  for al in results.I.split(","):
    ignoredAlleles[al] = 0

  
  # associate each protein (ENSP*)
  # with its associated GVP alleles
  d = defaultdict(list)

  # protein name -> chromosome
  prot2chrom = {}

  if results.D:
    dumpProts(protFile, csqFile)
    return 0
  
  with open(gvpPanelFile, 'r') as gvps:
    reader = csv.DictReader(gvps, delimiter="\t")
    for row in reader:

      peptide = uniquify(row[peptideCol])
      if peptide in ignoredAlleles:
        ignoredAlleles[peptide]=1 # mark that this allele was found in the gvp panel.
        continue
      
      start = uniquify(row[peptideStart])
      if type(start) != str:
        print("Parse error: the start coordinate is not singular!", start, file=sys.stderr, sep="\n")
        return 1

      # -1 converts from 1-based to 0/1 half open coordinates
      start = int(start)-1
    
      # when there's multiple SNPs in the peptide, the stop coordinates are a list separated by ;
      # it needs to just be one number though...
      stop = uniquify(row[peptideStop])
      if type(stop) != str:
        print("Parse error: the stop coordinate is not singular!", stop, file=sys.stderr, sep="\n")
        return 1


      stop = int(stop)
      # start/stop pair are in half-open coordinates
    
      prot = uniquify(row[protCol])
  

      prot2chrom[prot] =  row["chromosome"] 
      # This assumes a complete digest... lets make it more general
      #peptides = digest(peptide, asList=True)
      #if len(peptides) > 1:
        #print("Digest failed with: " , peptide, peptides, file=sys.stderr, sep="\n")
        #return 1
      
      al = Allele(prot, start, stop, isoleucine2leucine(peptide), peptide)
      d[prot].append(al)


  lut = open(results.L, "w")
  print("ProteinID", "Start", "Stop", "OriginalStart", "OriginalStop", "Peptide", sep="\t", file=lut)
  # housekeeping-- let's ensure that the allele calls we ask for are actually being recovered by the exome sequencing
  # if not, we may want to alter the GVP panel
  marginalAlleleCounts = dict()
  # changes d:
  # before, protein ID -> list of alleles
  # now protein ID -> list of loci (tuples of alleles)
  for prot, alleles in d.items():
    alleles.sort(key=lambda allele: (allele.protStart, allele.protStop))
    for a in alleles:
      marginalAlleleCounts[ a.seq ] = 0
    d[prot] = makeLoci(alleles)
    for loc in d[prot]:
      for allele in loc[1:]: # +1 to put it back into 1-based inclusive indexing (as per the input)
        print(prot, loc[0].protStart, loc[0].protStop, allele.protStart+1, allele.protStop, allele.seq, sep="\t", file=lut)

  lut.close()
  segsiteDict = {}

  # protein ID -> sequence
  # and
  # protein transcript ID -> protein ID
  # and
  # protein ID -> protein transcript ID
  (seqsDict, trans2prot, prot2trans) = getDictOfSeqs(protFile, d)

  # this is the PROTEIN haplotypes for
  # each distinct protein assayed
  # the protein haplotypes are difference encodings a la bcftools csq
  (haps, ids) = getHaplotypes(seqsDict, trans2prot, csqFile)

  # sort the peptides by their coordinate position (within the peptide)
  # this allows for detection of (potential) heterozygous sites...

  print("ID" , "Chromosome", "ProteinID", "ProteinStart", "ProteinStop" ,\
        "DetectableAllele1", "DetectableAllele2", "InferredAllele1", "InferredAllele2", \
        sep="\t")

  for prot, loci in d.items():
  
    trans = prot2trans[prot]
    protSeq = seqsDict[prot]
    transcript = prot2trans[prot]
    chromosome = prot2chrom[prot]
  
    if prot in haps:
      for locus in loci:
        
        refSeq = getRefAllele(locus[0], protSeq, True)
        peptides = getPeptideHaps( haps[prot], locus, refSeq, protSeq)
        refAllele = getConsistentAllele(locus, refSeq)
        refSeq = "".join(refSeq)
        for ID in ids:
          # nd stands for non-detect.
          # it's placed in lower case so that it is inconsistent with the amino acid alphabet
          alleles = ["nd", locus[0].seq, "nd", locus[0].seq]
          (b, e) = (locus[0].protStart+1, locus[0].protStop)
          # if a sample has no variants in the protein then they're not present
          # in the peptides data structure. As per the VCF file format,
          # this means that they're the reference sequence.
          if ID + "_1" in peptides:
            mAllele = peptides[  ID + "_1" ] # maternal allele 
          else:
            mAllele = (refSeq, refAllele)
          
          if mAllele[1].seq == 'null':
            alleles[0] = "nd"
          else:
            alleles[0] = mAllele[1].seq
          
          alleles[1] = mAllele[0]
            
          if ID + "_2" in peptides:
            pAllele = peptides[  ID + "_2" ] # and paternal allele
          else:
            pAllele = (refSeq, refAllele)
            
          if pAllele[1].seq == 'null':
            alleles[2] = "nd"
          else:
            alleles[2] = pAllele[1].seq
          
          alleles[3] = pAllele[0]

          print(ID, chromosome, prot, b, e,\
                alleles[0], alleles[2], alleles[1], alleles[3],\
                sep="\t")

          # let's evaluate the overall coverage of the peptide alleles
          if alleles[0] in marginalAlleleCounts:
            marginalAlleleCounts[alleles[0]] += 1
          if alleles[2] in marginalAlleleCounts:
            marginalAlleleCounts[alleles[2]] += 1

  # error checking
  marginalErrors=0
  for a,c in marginalAlleleCounts.items():
    if c == 0:
      print("Peptide allele:", "\t", a , \
            "was not found in the consequences file. Make sure this allele has sufficient coverage in the exome if you are to use it", \
            " otherwise the allele frequency will be taken as 5/2n in the RMP calculations (which may not be appropriate)", \
            file=sys.stderr, sep="\n");
      marginalErrors += 1

  for al,c in ignoredAlleles.items():
    if c==0: # this allele was not found (found==1) in the gvp panel
      print("Allele ", al , \
            "was asked to be dropped from the gvp panel... but it's not in the gvp panel. Not good!" ,\
            file=sys.stderr, sep="\n")
      marginalErrors += 1


      
  return marginalErrors

if __name__ == "__main__":
  sys.exit(main(sys.argv))

