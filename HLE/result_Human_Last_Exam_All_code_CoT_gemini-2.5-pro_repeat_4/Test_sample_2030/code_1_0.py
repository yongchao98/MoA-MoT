def get_reverse_complement(dna_sequence):
    """Calculates the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATGC", "TACG")
    return dna_sequence.upper().translate(complement_map)[::-1]

# Based on the 6-frame translation analysis, we identified Frame 6 as the target.
# The DNA for Frame 6 is a segment of the reverse complement of the original sequence.
# Reverse complement: 5'-ACCACTTGTGCGGGGAAG-3'
# The codons for Frame 6 are: CAC TTG TGC GGG GAA
# The amino acids are:          H   L   C   G   E

# The two amino acids fitting the SNP criteria are Leucine (L, non-polar) and Glutamic Acid (E, polar).
original_codons = ["CAC", "TTG", "TGC", "GGG", "GAA"]
original_aa = ["Histidine", "Leucine", "Cysteine", "Glycine", "Glutamic Acid"]

# SNP 1: A polar amino acid (Glutamic Acid, E) is changed to a stop codon.
# The codon for E is GAA. A single SNP (G->T) changes it to TAG (Stop).
modified_e_codon = "TAG"

# SNP 2: A non-polar amino acid (Leucine, L) is changed to a Cysteine.
# The codon for L is TTG. A single SNP (G->T) changes it to TGT (Cysteine).
modified_l_codon = "TGT"

# The modified sequence of codons is therefore:
modified_codons = ["CAC", modified_l_codon, "TGC", "GGG", modified_e_codon]

# The oligo binds to the sequence that is translated into amino acids.
# We exclude the stop codon (TAG) at the end.
target_dna_sequence = "".join(modified_codons[:-1])

# The oligo sequence is the reverse complement of this target sequence.
oligo_sequence = get_reverse_complement(target_dna_sequence)

print("1. Original codons in target frame (Frame 6):", " ".join(original_codons))
print("2. Corresponding amino acids:", " -> ".join(original_aa))
print("3. Modified codons after SNPs:", " ".join(modified_codons))
print("4. DNA sequence for oligo to bind (excluding stop codon):")
print(f"   5' {target_dna_sequence} 3'")
print("\n5. Final oligo sequence (reverse complement of the target):")
print(f"   5' {oligo_sequence} 3'")
