import textwrap

def get_reverse_complement(dna):
    """Computes the reverse complement of a DNA string."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, 'N') for base in reversed(dna))

def translate_dna(dna, frame):
    """Translates a DNA sequence into a peptide sequence."""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }
    
    peptide = []
    codons = []
    
    seq_in_frame = dna[frame-1:]
    
    for i in range(0, len(seq_in_frame), 3):
        codon = seq_in_frame[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
            peptide.append(codon_table.get(codon, 'X'))
            
    return "".join(peptide), codons

def main():
    original_seq = "CTTCCCCGCACAAGTGGT"
    rev_comp_seq = get_reverse_complement(original_seq)

    print("--- Step 1: Translating all 6 reading frames ---")
    
    frames = {}
    # Forward frames
    frames[1] = translate_dna(original_seq, 1)
    frames[2] = translate_dna(original_seq, 2)
    frames[3] = translate_dna(original_seq, 3)
    # Reverse complement frames
    frames[4] = translate_dna(rev_comp_seq, 1)
    frames[5] = translate_dna(rev_comp_seq, 2)
    frames[6] = translate_dna(rev_comp_seq, 3)
    
    print(f"Original Sequence: 5' {textwrap.fill(original_seq, width=18)} 3'")
    for i in range(1, 4):
        print(f"Frame +{i}: {' '*(i-1)}{frames[i][0]}")

    print(f"\nReverse Complement: 5' {textwrap.fill(rev_comp_seq, width=18)} 3'")
    for i in range(4, 7):
        print(f"Frame -{i-3}: {' '*(i-4)}{frames[i][0]}")

    print("\n--- Step 2: Finding the frame with two unique amino acids ---")

    all_aas = set("".join(f[0] for f in frames.values()))
    target_frame_id = -1
    unique_aas_in_target = []
    
    for i, (peptide, codons) in frames.items():
        other_peptides = "".join(frames[j][0] for j in frames if i != j)
        unique_aas = [aa for aa in set(peptide) if aa not in other_peptides]
        if len(unique_aas) == 2:
            target_frame_id = i
            unique_aas_in_target = unique_aas
            break
            
    print(f"Identified Target Frame: Frame {target_frame_id}")
    print(f"Peptide Sequence: {frames[target_frame_id][0]}")
    print(f"Unique Amino Acids in this frame: {unique_aas_in_target[0]} and {unique_aas_in_target[1]}")

    print("\n--- Step 3: Applying SNP mutations ---")
    
    # Per problem statement and analysis, F must be non-polar and W must be treated as polar.
    # Non-polar to Cysteine (Cys) mutation
    non_polar_aa = 'F'
    original_f_codon = "TTC"
    modified_c_codon = "TGC" # TTC -> TGC (T->G at pos 2)
    print(f"Mutation 1 (Non-polar to Cys): The non-polar amino acid Phenylalanine (F), coded by {original_f_codon}, is changed to Cysteine (C), coded by {modified_c_codon}.")

    # Polar to Stop mutation
    polar_aa = 'W'
    original_w_codon = "TGG"
    modified_stop_codon = "TGA" # TGG -> TGA (G->A at pos 3)
    print(f"Mutation 2 (Polar to Stop): The polar amino acid Tryptophan (W), coded by {original_w_codon}, is changed to a Stop codon (*), coded by {modified_stop_codon}.")

    # Reconstruct the modified DNA sequence for the frame
    original_codons = frames[target_frame_id][1]
    modified_codons = []
    for codon in original_codons:
        if codon == original_f_codon:
            modified_codons.append(modified_c_codon)
        elif codon == original_w_codon:
            modified_codons.append(modified_stop_codon)
        else:
            modified_codons.append(codon)

    print("\n--- Step 4: Designing the oligo ---")
    print(f"Original translated DNA sequence: 5'-{' '.join(original_codons)}-3'")
    print(f"Modified DNA sequence:          5'-{' '.join(modified_codons)}-3'")

    # The oligo targets the sequence coding for amino acids, which is up to the stop codon
    oligo_target_dna = "".join(modified_codons).split(modified_stop_codon)[0]
    
    print(f"\nThe oligo will bind to the sequence coding for amino acids: 5'-{oligo_target_dna}-3'")
    
    final_oligo_sequence = get_reverse_complement(oligo_target_dna)
    print(f"\nThe final oligo sequence is the reverse complement of this target.")
    print(f"Final Oligo Sequence (5' to 3'): {final_oligo_sequence}")
    
    return final_oligo_sequence

final_answer = main()
print(f"\n<<<{final_answer}>>>")