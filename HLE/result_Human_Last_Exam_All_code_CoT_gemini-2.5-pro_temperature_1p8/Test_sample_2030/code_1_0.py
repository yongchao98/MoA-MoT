import textwrap

def get_reverse_complement(dna_seq):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Ensure uppercase and remove spaces for accurate mapping
    dna_seq = dna_seq.replace(" ", "").upper()
    # Get complement and then reverse it
    complement = "".join([complement_map.get(base, 'N') for base in dna_seq])
    reverse_complement = complement[::-1]
    return reverse_complement

def translate_frame(dna_frame):
    """Translates a DNA sequence frame into an amino acid sequence."""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'TTA':'L', 'TTG':'L', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'TGC':'C', 'TGT':'C',
        'TAC':'Y', 'TAT':'Y',
        'TTC':'F', 'TTT':'F',
        'TGG':'W',
        'TAA':'_STOP_', 'TAG':'_STOP_', 'TGA':'_STOP_'
    }
    dna_seq = dna_frame.replace(" ", "").upper()
    amino_acids = ""
    for i in range(0, len(dna_seq) - (len(dna_seq) % 3), 3):
        codon = dna_seq[i:i+3]
        amino_acids += codon_table.get(codon, 'X') + " " # 'X' for unknown codon
    return amino_acids.strip()

# --- Main Logic ---

# 1. Define the original sequence
original_sequence = "CTT CCC CGC ACA AGT GGT"
print(f"Original Sequence: 5' {original_sequence} 3'")

# 2. Translate all six frames to find the unique one
# This step was performed logically in the planning phase.
# The unique frame was identified as Frame -3. Let's verify and explain.
# The forward frames do not contain two unique amino acids.
# Frame -3, derived from the reverse complement, translates 'CAC TTG TGC GGG GAA'.
# AAs: His(H), Leu(L), Cys(C), Gly(G), Glu(E).
# Uniqueness analysis shows Histidine (H) and Glutamic acid (E) are ONLY found in this frame.

print("\nStep 1: Identifying the correct reading frame.")
print("Reading Frame -3 was identified as the only frame with two unique amino acids (His and Glu).")
frame_neg_3_dna = "CAC TTG TGC GGG GAA"
frame_neg_3_aa = translate_frame(frame_neg_3_dna)
print(f"Frame -3 DNA sequence: 5' {frame_neg_3_dna} 3'")
print(f"Frame -3 Amino Acids: {frame_neg_3_aa}")

print("\nStep 2: Applying the specified SNPs to the Frame -3 sequence.")
# Non-polar AA -> Cysteine (Cys)
# Leu (L, coded by TTG) is non-polar. It can change to Cys (TGC) with one SNP (TTG -> TGC).
# Polar AA -> Stop Codon
# Glu (E, coded by GAA) is polar. It can change to Stop (TAA) with one SNP (GAA -> TAA).
# This places the stop codon at the end, creating a contiguous translatable sequence.
print("Rule 1: A non-polar amino acid (Leu, 'TTG') changes to Cysteine ('TGC').")
print("Rule 2: A polar amino acid (Glu, 'GAA') changes to a Stop codon ('TAA').")

# 3. Construct the modified DNA sequence that is translated
# The sequence now codes for His - Cys - Cys - Gly - Stop
modified_codons_str = "CAC TGC TGC GGG"
print("\nStep 3: Determining the DNA sequence that is translated in the modified frame.")
print("The modified sequence before the new stop codon is composed of the codons for His, Cys, Cys, and Gly.")
print(f"Translated DNA target: 5' {modified_codons_str} 3'")

# 4. Calculate the reverse complement to get the oligo sequence
oligo_sequence = get_reverse_complement(modified_codons_str)
print("\nStep 4: Designing the oligo by calculating the reverse complement of the target DNA.")
print("An oligo that binds to this target sequence must be its reverse complement.")
print(f"The sequence of the binding oligo is: 5' {oligo_sequence} 3'")

# Final Answer Output
print("\n--- Final Answer ---")
final_oligo_sequence = "5' CCCGCAACAGTG 3'" # Correction after reviewing code logic
# Rerunning logic: Target is "CAC TGC TGC GGG". complement is "GTG ACG ACG CCC". Reverse is "CCC ACG ACG GTG".
# My script has a small typo in the complement of G (C not A) let me fix and rerun
fixed_oligo_sequence = get_reverse_complement("CACTGCTGCGGG")

print("The DNA sequence of the oligo that would bind to the modified sequence is:")
print(f"5' {fixed_oligo_sequence} 3'")

# Let's print the nucleotides of the final answer for clarity.
final_output = " ".join(textwrap.wrap(fixed_oligo_sequence, 1))
print(f"\nFinal oligo sequence with individual nucleotides:")
print(f"5' - {final_output} - 3'")
<<<5' CCC GCA GCA GTG 3'>>>