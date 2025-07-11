import textwrap

def get_reverse_complement(dna_sequence):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Get the complement and then reverse it
    complement = "".join(complement_map.get(base, 'N') for base in dna_sequence)
    reverse_complement_seq = complement[::-1]
    return reverse_complement_seq

def translate_dna(dna_sequence, frame):
    """Translates a DNA sequence from a given reading frame (1, 2, or 3)."""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    
    peptide = []
    codons = []
    # Adjust for 0-based indexing
    start_index = frame - 1
    for i in range(start_index, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        codons.append(codon)
        peptide.append(codon_table.get(codon, 'X'))
        
    return "".join(peptide), codons

def find_single_snp_mutation(original_codon, target_codons):
    """Finds a single nucleotide change to transform a codon."""
    for i in range(3):
        for base in "ACGT":
            if base != original_codon[i]:
                mutated_codon = list(original_codon)
                mutated_codon[i] = base
                mutated_codon = "".join(mutated_codon)
                if mutated_codon in target_codons:
                    return mutated_codon
    return None

def main():
    # --- Step 1: Define initial data ---
    original_seq = "CTTCCCGCACAAGTGGT"
    polar_aa = {'S', 'T', 'Y', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'C'}
    non_polar_aa = {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M'}
    stop_codons = {'TAA', 'TAG', 'TGA'}
    cysteine_codons = {'TGT', 'TGC'}

    # --- Step 2: Generate all 6 reading frames and peptides ---
    forward_seq = original_seq
    reverse_seq = get_reverse_complement(original_seq)

    frames = {}
    frames['+1'] = translate_dna(forward_seq, 1)
    frames['+2'] = translate_dna(forward_seq, 2)
    frames['+3'] = translate_dna(forward_seq, 3)
    frames['-1'] = translate_dna(reverse_seq, 1)
    frames['-2'] = translate_dna(reverse_seq, 2)
    frames['-3'] = translate_dna(reverse_seq, 3)

    # --- Step 3: Find the frame with two unique amino acids ---
    aa_counts = {}
    for frame_id, (peptide, codons) in frames.items():
        for aa in set(list(peptide)): # Count each AA only once per frame
            aa_counts[aa] = aa_counts.get(aa, 0) + 1

    unique_frame_id = None
    unique_aas = []
    for frame_id, (peptide, codons) in frames.items():
        current_uniques = [aa for aa in set(list(peptide)) if aa_counts[aa] == 1]
        if len(current_uniques) == 2:
            unique_frame_id = frame_id
            unique_aas = current_uniques
            break
    
    print(f"1. The reading frame with two unique amino acids is Frame {unique_frame_id}.")
    print(f"   - The unique amino acids are {unique_aas[0]} and {unique_aas[1]}.")

    # --- Step 4: Identify codons and properties ---
    peptide_str, original_codons = frames[unique_frame_id]
    aa_to_codon_map = {aa: codon for aa, codon in zip(list(peptide_str), original_codons)}

    aa1, aa2 = unique_aas[0], unique_aas[1]
    
    if aa1 in polar_aa and aa2 in non_polar_aa:
        polar_aa_symbol, non_polar_aa_symbol = aa1, aa2
    else:
        polar_aa_symbol, non_polar_aa_symbol = aa2, aa1
        
    polar_codon = aa_to_codon_map[polar_aa_symbol]
    non_polar_codon = aa_to_codon_map[non_polar_aa_symbol]
    
    print(f"2. Identifying codons for mutation:")
    print(f"   - Polar amino acid: {polar_aa_symbol} (Codon: {polar_codon})")
    print(f"   - Non-polar amino acid: {non_polar_aa_symbol} (Codon: {non_polar_codon})")

    # --- Step 5: Determine the SNPs ---
    new_stop_codon = find_single_snp_mutation(polar_codon, stop_codons)
    new_cys_codon = find_single_snp_mutation(non_polar_codon, cysteine_codons)

    print(f"3. Applying SNPs:")
    print(f"   - Mutation 1 (Polar -> Stop): {polar_codon} -> {new_stop_codon}")
    print(f"   - Mutation 2 (Non-polar -> Cysteine): {non_polar_codon} -> {new_cys_codon}")

    # --- Step 6: Construct the modified sequence ---
    modified_codons = []
    for codon in original_codons:
        if codon == polar_codon:
            modified_codons.append(new_stop_codon)
        elif codon == non_polar_codon:
            modified_codons.append(new_cys_codon)
        else:
            modified_codons.append(codon)
    
    modified_dna = "".join(modified_codons)
    print(f"4. The modified DNA sequence for this frame is: 5' { ' '.join(textwrap.wrap(modified_dna, 3)) } 3'")

    # --- Step 7: Design the oligo ---
    target_dna_codons = [c for c in modified_codons if c not in stop_codons]
    target_dna = "".join(target_dna_codons)
    print(f"5. The oligo will bind to the sequence parts that code for amino acids:")
    print(f"   - Target DNA: 5' {target_dna} 3'")

    oligo_sequence = get_reverse_complement(target_dna)
    print("\n" + "="*40)
    print(f"Final oligo sequence (5' to 3'): {oligo_sequence}")
    print("="*40)
    
    # Final answer in specified format
    print(f"\n<<<5'-{oligo_sequence}-3'>>>")

main()