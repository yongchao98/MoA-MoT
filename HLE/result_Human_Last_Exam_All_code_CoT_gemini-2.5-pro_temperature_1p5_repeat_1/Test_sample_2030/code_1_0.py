import collections

def get_reverse_complement(dna_sequence):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, 'N') for base in reversed(dna_sequence))

def translate_dna(dna_sequence):
    """Translates a DNA sequence into an amino acid sequence."""
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    amino_acids = []
    for i in range(0, len(dna_sequence) - (len(dna_sequence) % 3), 3):
        codon = dna_sequence[i:i+3]
        amino_acids.append(codon_table.get(codon, '?'))
    return amino_acids

def solve_oligo_puzzle():
    """Main function to solve the bioinformatics puzzle."""
    
    # Step 1: Define initial data
    original_seq = "CTTCCCGCACAAGTGGT"
    print(f"Original Sequence: 5' {original_seq} 3'\n")
    
    rev_comp_seq = get_reverse_complement(original_seq)

    # Step 2: Generate and translate all 6 reading frames
    frames_dna = {
        'F1': [original_seq[i:i+3] for i in range(0, len(original_seq), 3) if i+3 <= len(original_seq)],
        'F2': [original_seq[i:i+3] for i in range(1, len(original_seq), 3) if i+3 <= len(original_seq)],
        'F3': [original_seq[i:i+3] for i in range(2, len(original_seq), 3) if i+3 <= len(original_seq)],
        'F4': [rev_comp_seq[i:i+3] for i in range(0, len(rev_comp_seq), 3) if i+3 <= len(rev_comp_seq)],
        'F5': [rev_comp_seq[i:i+3] for i in range(1, len(rev_comp_seq), 3) if i+3 <= len(rev_comp_seq)],
        'F6': [rev_comp_seq[i:i+3] for i in range(2, len(rev_comp_seq), 3) if i+3 <= len(rev_comp_seq)],
    }
    
    frames_aa = {name: translate_dna("".join(codons)) for name, codons in frames_dna.items()}

    print("--- Reading Frames and Amino Acid Sequences ---")
    for name, codons in frames_dna.items():
        print(f"Frame {name}: {' '.join(codons)} -> {' '.join(frames_aa[name])}")
    
    # Step 3: Find the frame with two unique amino acids
    all_aa = [aa for seq in frames_aa.values() for aa in seq if aa != '_']
    aa_counts = collections.Counter(all_aa)
    
    unique_aa_per_frame = {}
    for name, seq in frames_aa.items():
        uniques = {aa for aa in seq if aa_counts[aa] == 1 and aa != '_'}
        if uniques:
            unique_aa_per_frame[name] = uniques

    special_frame_name = None
    for name, uniques in unique_aa_per_frame.items():
        if len(uniques) == 2:
            special_frame_name = name
            break
            
    print("\n--- Identifying the Special Frame ---")
    print(f"Found frame with two unique amino acids: {special_frame_name}")
    print(f"Unique amino acids in this frame are: {', '.join(unique_aa_per_frame[special_frame_name])}")

    # Step 4: Analyze the contradiction and state assumption
    # The unique AAs are F and W, both non-polar. The prompt states one changed AA was polar.
    # This is a contradiction. We proceed by assuming a typo and that a NON-POLAR AA was changed to STOP.
    print("\n--- Applying SNPs based on problem description (with one assumption) ---")
    print("Assumption: The problem contains a typo. It should state that a non-polar amino acid (not polar) was changed to a stop codon, as the unique amino acids (F and W) are both non-polar.")
    
    # Step 5: Determine the modifications
    # The SNPs must affect the codons for the unique amino acids F (TTC) and W (TGG)
    # Change 1: Non-polar (W, codon TGG) -> Stop codon. A single SNP (G->A) changes TGG to TGA (STOP).
    # Change 2: Non-polar (F, codon TTC) -> Cysteine. A single SNP (T->G) changes TTC to TGC (Cysteine).
    original_orf_codons = frames_dna[special_frame_name]
    modified_orf_codons = []
    for codon in original_orf_codons:
        if codon == "TTC": # Original F codon
            modified_orf_codons.append("TGC") # Becomes Cys
        elif codon == "TGG": # Original W codon
            modified_orf_codons.append("TGA") # Becomes STOP
        else:
            modified_orf_codons.append(codon)
    
    print("\nOriginal Frame 2 DNA: ", " ".join(original_orf_codons))
    print("Modified Frame 2 DNA: ", " ".join(modified_orf_codons))

    # Step 6: Identify the translated sequence and design the oligo
    translated_part = []
    for codon in modified_orf_codons:
        if translate_dna(codon)[0] == '_':
            break
        translated_part.append(codon)
    
    target_dna = "".join(translated_part)
    oligo_sequence = get_reverse_complement(target_dna)
    
    print("\n--- Designing the Oligo ---")
    print(f"The DNA sequence to be translated in the modified frame is: 5' {target_dna} 3'")
    print(f"The oligo that binds to this sequence (its reverse complement) is:")
    print(f"Final Oligo Sequence: 5' {oligo_sequence} 3'")
    
    return oligo_sequence

final_answer = solve_oligo_puzzle()
print(f"\n<<<5'-{final_answer}-3'>>>")