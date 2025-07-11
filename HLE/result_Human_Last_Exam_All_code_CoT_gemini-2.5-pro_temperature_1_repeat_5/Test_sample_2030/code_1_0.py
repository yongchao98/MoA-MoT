def get_reverse_complement(dna_sequence):
    """
    Calculates the reverse complement of a DNA sequence.
    """
    complement_map = str.maketrans('ATCG', 'TAGC')
    # Complement the sequence, then reverse it
    return dna_sequence.upper().translate(complement_map)[::-1]

# Step 1: Define the original DNA sequence and its reverse complement.
original_sequence_5_to_3 = "CTTCCCGCACAAGTGGT"
# The reverse complement is the sequence read backwards and complemented.
# 5'-CTTCCCGCACAAGTGGT-3' -> Reverse: 3'-TGGTGACACGCCC TTC-5' -> RC: 5'-ACCACTTGTCGCGGGAAG-3'
# This leads to Frame -1, which is the frame of interest.
frame_minus_1_dna = "ACCACTTGTCGCGGGAAG"

# Step 2 & 3: Identify the codons for the unique amino acids based on prior analysis.
# Frame -1 translation: T T C R G K
# Unique amino acids are C (Cysteine) and K (Lysine).
# Codon for Cysteine (C) in this frame is TGT.
# Codon for Lysine (K) in this frame is AAG.
# We assume Cysteine is non-polar and Lysine is polar to solve the problem.

# Step 4: Determine the modified sequence based on the SNP rules.
# The polar AA (Lysine, AAG) changes to a stop codon. A single SNP A->T gives TAG.
# The non-polar AA (Cysteine, TGT) changes to Cysteine. A single SNP T->C gives TGC.
original_lysine_codon = "AAG"
modified_lysine_codon = "TAG" # Becomes a stop codon

original_cysteine_codon = "TGT"
modified_cysteine_codon = "TGC" # Synonymous change to Cysteine

# The original Frame -1 sequence was composed of these codons:
# ACC ACT TGT CGC GGG AAG
# We replace TGT with TGC and AAG with TAG to get the modified sequence.
modified_frame_dna = "ACCACTTGCCGCGGGTAG"

# Step 5: Design the oligo.
# The oligo should bind to the part of the modified sequence that is translated.
# The stop codon (TAG) is not translated into an amino acid.
oligo_target_sequence = modified_frame_dna[:-3] # Remove the last 3 bases (the stop codon)

# The oligo sequence is the reverse complement of its target.
final_oligo_sequence = get_reverse_complement(oligo_target_sequence)

print("The DNA sequence of the oligo is:")
print(f"5' {final_oligo_sequence} 3'")

# The final answer format requires printing the sequence itself.
# Let's print it in the requested format.
# We will print the final oligo sequence base by base as requested.
final_equation = ""
for base in final_oligo_sequence:
    final_equation += base + " "
print("Final Equation: 5' " + final_equation.strip() + " 3'")
# The final result is the sequence itself.
#<<<5' CCC GCG GCA AGT GGT 3'>>>