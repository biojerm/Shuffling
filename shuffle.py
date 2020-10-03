#!/usr/bin/env python

import sys
import re
import difflib

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import pairwise2  # allows us to do pairwise alignment on our proteins
from Bio.Seq import Seq


def align_with_blosum62(aa_seq1, aa_seq2):
    """
    Creates a protein alignment when given two amino acid sequences.  Tuple of top alignments is returned
    """

    # note: depending on the sequence homology it may make sense to use
    # another blosum matrix (or different gap_open gap_close)
    # matrix = matlist.align_with_blosum62
    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -12  # cost to open a gap
    gap_extend = -3  # cost to extend a gap

    alignments = pairwise2.align.globalds(
        aa_seq1, aa_seq2, matrix, gap_open, gap_extend)

    return(alignments[0])


def find_alignment_gaps(alignment):
    """
    Returns location (in list form) of alignment gaps in a sequence
    """
    return([char.start() for char in re.finditer("-", alignment)])


def insert_char(mystring, position, chartoinsert):
    mystring = mystring[:position] + chartoinsert + mystring[position:]
    return mystring


def insert_nt_gaps(dna_sequence, alignment_gaps):
    for position in alignment_gaps:
        dna_sequence = insert_char(dna_sequence, position * 3, '---')
    return(dna_sequence)


def calculate_mutation_spaces(mutation_list):
    mutation_list.insert(0, 0)  # add zero to front of list
    mutation_location = []
    for e, m in enumerate(mutation_list):
        mutation_location.append(m - mutation_list[e - 1])
        mutation_location[0] = 0  # I add in the zero here to make sure the first mutation in the list is counted (its baggage is removed by slicing the results from left_right_greater)
    # print(mutation_location)
    return(mutation_location)
# calculate_mutation_spaces(dist)


def left_right_greater(mutation_space_distances):
    """
    based on distances from mutation to mutation, this will return the direction with more homology.
    Will return either left or right (probably should be upstream/downstream)
    """
    direction_of_greater_homology = []
    for e, distance in enumerate(mutation_space_distances):
        if e == 0:
            direction_of_greater_homology.append('right')
        elif e == len(mutation_space_distances) - 1:  # I think i need to have -1 b/c index start at 0 double check I am not forgetting the last mutation
            direction_of_greater_homology.append('left')
        elif distance > mutation_space_distances[e + 1]:
            direction_of_greater_homology.append('left')
        elif distance < mutation_space_distances[e + 1]:
            direction_of_greater_homology.append('right')
        else:  # if homologies are equal just pick left because...choose one.
            direction_of_greater_homology.append('left')
    return(direction_of_greater_homology[1:])  # see note in calculate_mutation_spaces for explanation of slice


acid_dic = {
    'M': ['ATG'],
    'W': ['TGG'],
    'F': ['TTT', 'TTC'],
    'Y': ['TAT', 'TAC'],
    'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['TGT', 'TGC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']}


def similar_codons(ref_codon, test_codon):
    test_aa = str(Seq(test_codon).translate())
    difference_ratio = []
    tested_codons = []
    for codon in acid_dic[test_aa]:
        dif = difflib.SequenceMatcher(None, ref_codon.upper(), codon.upper())
        difference_ratio.append(dif.ratio())
        tested_codons.append(codon)
    max_homology = max(difference_ratio)
    best_matches = []
    for e, ratio in enumerate(difference_ratio):
        if ratio == max_homology:
            best_matches.append(tested_codons[e])
    if test_codon.upper() in best_matches:   # If the codon max homology, there is no reason to change it
        best_matches = [test_codon.upper()]
# Write in something about how if original codon is same homology than other codon
# return original codon.  No reason to change what is essentially equivalent.
    return(best_matches)


def codon_analyzer(ref_codon, codon_list):
    similarity_dict = {}
    for codon in codon_list:
            similarity_score = []
            for e, v in enumerate(codon):
                if v.upper() == ref_codon[e].upper():
                    similarity_score.append('1')
                else:
                    similarity_score.append('0')
            similarity_dict[codon] = int(''.join(similarity_score))  # combine the '0' and '1' then turn into a number
    return(similarity_dict)


def codon_switch(sequence, codon, position):
    first = sequence[:position + 1]
    last = sequence[position + 4:]
    insert = codon
    return(first + insert + last)

# ######  Start of the actual coding being run.  Functions go above this point


if __name__ == '__main__':
    fa_seq_records = []

    # Import files from given name in terminal window, alphabet is also assigned (not currently that useful of a feature)
    with open(sys.argv[1]) as file:
        for seq in SeqIO.parse(file, "fasta"):
            fa_seq_records.append(seq)
        file.close()

    # First entry in fasta file is the reference for codon harmonization
    sequence_reference = fa_seq_records[0]
    # the rest of the entries in the fasta file
    sequences_to_harmonize = fa_seq_records[1:]

    for record in sequences_to_harmonize:

        # Working DNA variables
        reference_dna = sequence_reference.seq
        dna_to_harmonize = record.seq

        # Working amino acid variables
        reference_aa = fa_seq_records[0].seq.translate(to_stop=True)
        aa_to_harmonize = record.seq.translate(to_stop=True)

        # alignment of amino acid sequences and also assignment of alignment score variables
        top_aln = align_with_blosum62(reference_aa, aa_to_harmonize)
        aln_ref, aln_harmonize, score, begin, end = top_aln

        # Inserting gaps in sequences based on alignments, so sequences are same length. Gaps are '---'.
        ref_seq_with_gaps = insert_nt_gaps(reference_dna, find_alignment_gaps(aln_ref))
        dna_to_harmonize_with_gaps = insert_nt_gaps(dna_to_harmonize, find_alignment_gaps(aln_harmonize))

        # Copies codons from reference sequence if amino acids match.
        # If gap is present in reference sequence, the codon is copied from the sequence being harmonized
        # lastly if the amino acid is different the codon from the haromonization sequence is copied and
        #    the location noted in mismatch_aa_locations
        new_harmonized_sequence = ''
        mismatch_aa_locations = []
        for i, aa in enumerate(aln_ref):
            if aa == aln_harmonize[i]:
                new_harmonized_sequence += ref_seq_with_gaps[(i * 3):(i * 3) + 3]
            elif aa == '-':
                new_harmonized_sequence += dna_to_harmonize_with_gaps[(i * 3):(i * 3) + 3]
            else:
                new_harmonized_sequence += dna_to_harmonize_with_gaps[(i * 3):(i * 3) + 3]
                mismatch_aa_locations.append(i)

        # add back stop codon (left off because of alignment).
        new_harmonized_sequence += dna_to_harmonize_with_gaps[-3:]

        """
        This is the part of the code I am least happy with. I am trying to optimize the codons for partial amino acids that are not the same
        using codon variations to get and extra NT or two more homology.

        The logic is to calculate the difference between two mutations and improve the homology on the side with longer similarity stretches
        """

        # Calculate if larger homology stretch is on the left or right of a sequence
        # IMPORTANT calculate_mutation_spaces adds a 0 to the mutation list to allow calculate the distance to the start codon
        #     as a result the first item of the list sliced out.
        left_right_list = left_right_greater(calculate_mutation_spaces(mismatch_aa_locations))

        for e, aa_position in enumerate(mismatch_aa_locations[1:]):
            nt_position = aa_position * 3
            ref = str(ref_seq_with_gaps[nt_position:nt_position + 3])  # these are being convereted to biopython sequences at some point to I need to turn them into strings
            target = str(new_harmonized_sequence[nt_position:nt_position + 3])
            codon_to_switch = ''
            if target == '---':
                pass
            elif len(set(codon_analyzer(ref, similar_codons(ref, target)).values())) == 1:
                # print('original codon')
                codon_to_switch = target
            elif left_right_list[e] == 'left':
                # print('max value')
                codon_to_switch = max(codon_analyzer(ref, similar_codons(ref, target)), key=lambda key: codon_analyzer(ref, similar_codons(ref, target))[key])
            elif left_right_list[e] == 'right':
                # print('min value')
                codon_to_switch = min(codon_analyzer(ref, similar_codons(ref, target)), key=lambda key: codon_analyzer(ref, similar_codons(ref, target))[key])

            else:
                codon_to_switch = target

            codon_switch(new_harmonized_sequence, codon_to_switch, nt_position)

        new_harmonized_sequence = str(new_harmonized_sequence).replace('-', '')
        print(record.description + '--Harmonized to-->' + sequence_reference.description)
        print(new_harmonized_sequence)
