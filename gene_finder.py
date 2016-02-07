# -*- coding: utf-8 -*-
"""
Gene Finder

@author: Tehya Stockman

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """

    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """

    #the following are the only combinations of nucleotides allowed

    if nucleotide == 'A':
    	return 'T'
    elif nucleotide == 'G':
    	return 'C'
    elif nucleotide == 'T':
    	return 'A'
    elif nucleotide == 'C':
    	return 'G'
    

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("TATGCCC")  #This doc test is good, as as long as each of the characters in the string is giving the complement and in reverse order, the function is working properly
    'GGGCATA'
    """
   
    reverse = dna[::-1]                  #reverses the DNA strand
    
    new_reverse = ''                     #set new strand to empty set
    for char in reverse:                 #for each character in the strand
    	t = get_complement(char)         #uses function above to switch nucleotide
    	new_reverse = new_reverse + t    #adds current onto previous
  
    	
    return new_reverse
    
def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATGAG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGTAACTGAGA")  #Should not cut off at things that look like stop codons that are in the wrong frame, needs to be in sets of 3, make sure it won't malfunction past the first time it reaches what could be a stop codon, but will reccognize later ones
    'ATGAGTAAC'
    """
    # stop codons include: TAG, TAA, TGA

    temp = ''   #empty string for dna
   
    a = 'TAA'   #the three stop codons
    b = 'TAG'
    c = 'TGA'

    for i in range(0, len(dna), 3):
    	nucleotide = dna[i:i+3]                                    #range of each dna nucleotide
    	if nucleotide == a or nucleotide == b or nucleotide == c:  #each or must be followed by a complete statement
    		return temp                                            #returns the temporary list of nucleotides
    	temp += nucleotide                                         #updating list after makes it not add last nucleotide

    return temp
    

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        #The second test should make sure that function is only finding the ORFs of a single frame. 
        The strand is similar to the one above, but with an extra letter after the stop codon so 
        that the second strand of the first test should not show in the second.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGAATGTGCCC")
    ['ATGCATGAATGTAGA']
    """
    list_ORFs = []   #empty list for ORFs to be filled
    i = 0            #index to search through the list
    start_codon = 'ATG'

    while i < len(dna):     #beginning to search through list
    	nucleotide = dna[i:i+3]                      #index for every 3 nucleotides
    	if nucleotide == start_codon:
    		rest_of_strand = rest_of_ORF(dna[i:])   #goes from when find start codon to the end of the strand, calls function to get rid of stop codon
    		list_ORFs.append(rest_of_strand)        #adds to the list of ORFs

    	i = i+3  #updates the index by three, as amino acides can only be found every three nucleotides
   
    return list_ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        The doc test below is good because there are many ATG's and many of those are 
        nested within eachother. It covers the bases.

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []  #open a list for all the ORFs
    frame_1 = find_all_ORFs_oneframe(dna)    #ORFs for frame 1
    frame_2 = find_all_ORFs_oneframe(dna[1:])  #'' for frame 2
    frame_3 = find_all_ORFs_oneframe(dna[2:])  #'' for frame 3
                                 #concatenating all of the lists into a single list
    all_ORFs.extend(frame_1)
    all_ORFs.extend(frame_2)
    all_ORFs.extend(frame_3)
   

    return all_ORFs


    #find everything that begins with ATG, have to look at the multiples of three but adjust the frames


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    all_ORFs = []        #opens a list to put the strands

    other_strand = get_reverse_complement(dna)  #get reverse of the dna strand to find ORFs on other side
    strand_1 = find_all_ORFs(dna)
    strand_2 = find_all_ORFs(other_strand)
                                               #concatenating the lists into one
    all_ORFs.extend(strand_1)
    all_ORFs.extend(strand_2)

    return all_ORFs
    

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        This doc test works. I also tested it by printing the strands list to see
        all the different initial strands of the list from the previous function.
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    strands = find_all_ORFs_both_strands(dna)
    if strands == []:
    	return ''
    return max(strands, key = len)  #finds maximum strand based on length
    

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 

        I tested this function by printing the shuffled strings, as well as the maximum 
        string length, as both should vary each time the function is called.
    """
    
    list_sequences = []

    #find the longest ORF of shuffled dna and concatenate a list of them
    for i in range(num_trials):             
        x = longest_ORF(shuffle_string(dna))
    	list_sequences.append(x)

    #sort the list by length and find the longest strand of the list	
    list_sequences.sort(key = len)
    list_sequences
    return len(list_sequences[-1])
    
 
def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        These doc tests are sufficient, because in order to pass the second test,
        you must break up the dna into chunks of three, for which it will work.

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino_acids = ''
    for i in range(0, len(dna), 3):   #amino acids are made up of 3 nucleotides each
    	nucleotide = dna[i:i+3]  
        amino_acid = aa_table[nucleotide]   #calling the dictionary of amino acids to assign 
        amino_acids += amino_acid

    return amino_acids


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        Doc tests may not work for this, but one way to test that things are functioning
        is that all of the strings should start with M.
    """
    threshold = longest_ORF_noncoding(dna, 1500)  #randomly shuffled dna longest strand as threshold for possible protein length
    coding_ORF = find_all_ORFs_both_strands(dna) #list of dna strands
    list_protein_chain = []

    for i in coding_ORF:       #comparing non-shuffled dna to shuffled dna
        if len(i) > threshold:
    	    protein = coding_strand_to_AA(i)
    	    list_protein_chain.append(protein)

    print list_protein_chain

   
#if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose = True)
#get_reverse_complement('AATAC')
#rest_of_ORF('ATGCATGAATGTAGATAGATGTGCCC')
#find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
#find_all_ORFs("ATGCATGAATGTAG")
#find_all_ORFs_both_strands("ATGCATGAATGTAGATAGATGTGCCC")
#longest_ORF("ATGCATGAATGTAGATAGATGTGCCC")

from load import load_seq
dna = load_seq("./data/X73525.fa")
gene_finder(dna)
