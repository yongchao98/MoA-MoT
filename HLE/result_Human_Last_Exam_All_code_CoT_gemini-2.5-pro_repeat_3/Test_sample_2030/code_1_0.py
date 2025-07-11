import textwrap

def solve_oligo_puzzle():
    """
    Solves the bioinformatics puzzle by translating DNA, finding a unique reading frame,
    simulating SNPs, and designing a complementary oligo.
    """
    # Step 0: Define constants and helper functions
    original_seq = "CTT CCC CGC ACA AGT GGT".replace(" ", "")
    
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 
        'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 
        'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 
        'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 
        'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 
        'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    polar_aa = {'S', 'T', 'C', 'Y', 'N', 'Q', 'K', 'R', 'H', 'D', 'E'}
    non_polar_aa = {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M'}
    cysteine_codons = {'TGT', 'TGC'}
    stop_codons = {'TAA', 'TAG', 'TGA'}

    def reverse_complement(seq):
        complement_map = str.maketrans("ATCG", "TAGC")
        return seq.upper().translate(complement_map)[::-1]

    def translate_frame(seq):
        codons = textwrap.wrap(seq, 3)
        # Only consider full codons
        codons = [c for c in codons if len(c) == 3]
        amino_acids = [codon_table[c] for c in codons]
        return amino_acids, codons

    print(f"Original Sequence: 5' {original_seq} 3'\n")
    
    # Step 1: Translate all six reading frames
    print("--- Step 1: Translating all six reading frames ---")
    rev_comp_seq = reverse_complement(original_seq)
    print(f"Reverse Complement: 5' {rev_comp_seq} 3'")
    
    translations = {}
    codons_map = {}
    all_aa = set()
    
    # Forward frames
    for i in range(3):
        aa_seq, codon_seq = translate_frame(original_seq[i:])
        translations[f'Forward Frame {i+1}'] = aa_seq
        codons_map[f'Forward Frame {i+1}'] = codon_seq
        all_aa.update(aa_seq)

    # Reverse frames
    for i in range(3):
        aa_seq, codon_seq = translate_frame(rev_comp_seq[i:])
        translations[f'Reverse Frame {i-3}'] = aa_seq
        codons_map[f'Reverse Frame {i-3}'] = codon_seq
        all_aa.update(aa_seq)
    
    # Rename reverse frames for clarity (4, 5, 6)
    for i in range(3):
        key = f'Reverse Frame {i-3}'
        new_key = f'Reverse Frame {i+4}'
        translations[new_key] = translations.pop(key)
        codons_map[new_key] = codons_map.pop(key)

    for frame, aa_seq in translations.items():
        print(f"{frame}: {' '.join(aa_seq)}")
    print("-" * 20)

    # Step 2: Identify the frame with two unique amino acids
    print("\n--- Step 2: Identifying the unique frame ---")
    unique_frame = None
    unique_aa_info = []

    for frame, aa_seq in translations.items():
        other_aa = set()
        for other_frame, other_aa_seq in translations.items():
            if frame != other_frame:
                other_aa.update(other_aa_seq)
        
        current_unique_aa = []
        for i, aa in enumerate(aa_seq):
            if aa not in other_aa:
                current_unique_aa.append({'aa': aa, 'codon': codons_map[frame][i]})
        
        if len(current_unique_aa) == 2:
            unique_frame = frame
            unique_aa_info = current_unique_aa
            break
            
    print(f"Found unique frame: {unique_frame}")
    print(f"Unique Amino Acids and their codons: {unique_aa_info[0]['aa']}({unique_aa_info[0]['codon']}), {unique_aa_info[1]['aa']}({unique_aa_info[1]['codon']})")
    print("-" * 20)
    
    # Step 3: Verify amino acid properties and identify roles
    print("\n--- Step 3: Verifying amino acid properties ---")
    polar_candidate = None
    non_polar_candidate = None
    
    for info in unique_aa_info:
        if info['aa'] in polar_aa:
            polar_candidate = info
        elif info['aa'] in non_polar_aa:
            non_polar_candidate = info
            
    if not (polar_candidate and non_polar_candidate):
        print("Error: Could not find one polar and one non-polar unique amino acid.")
        return

    print(f"Polar amino acid to become STOP: {polar_candidate['aa']} (codon {polar_candidate['codon']})")
    print(f"Non-polar amino acid to become Cysteine: {non_polar_candidate['aa']} (codon {non_polar_candidate['codon']})")
    print("-" * 20)

    # Step 4: Simulate the SNPs
    print("\n--- Step 4: Simulating the SNPs ---")
    # Non-polar to Cysteine
    modified_non_polar_codon = None
    for cys_codon in cysteine_codons:
        diff = sum(1 for a, b in zip(non_polar_candidate['codon'], cys_codon) if a != b)
        if diff == 1:
            modified_non_polar_codon = cys_codon
            break
    print(f"SNP 1: {non_polar_candidate['codon']} ({non_polar_candidate['aa']}) -> {modified_non_polar_codon} (C)")

    # Polar to Stop
    modified_polar_codon = None
    for stop_codon in stop_codons:
        diff = sum(1 for a, b in zip(polar_candidate['codon'], stop_codon) if a != b)
        if diff == 1:
            modified_polar_codon = stop_codon
            break
    print(f"SNP 2: {polar_candidate['codon']} ({polar_candidate['aa']}) -> {modified_polar_codon} (*)")
    print("-" * 20)

    # Step 5: Construct the modified sequence and design the oligo
    print("\n--- Step 5: Designing the oligo ---")
    
    # Reconstruct the original codon sequence for the unique frame
    original_codons = codons_map[unique_frame]
    modified_codons = []
    
    # Replace the original codons with the modified ones
    for codon in original_codons:
        if codon == non_polar_candidate['codon']:
            modified_codons.append(modified_non_polar_codon)
        elif codon == polar_candidate['codon']:
            modified_codons.append(modified_polar_codon)
        else:
            modified_codons.append(codon)
            
    # Find the part that is translated into amino acids (before the stop codon)
    target_codons = []
    for codon in modified_codons:
        if codon in stop_codons:
            break
        target_codons.append(codon)

    target_sequence = "".join(target_codons)
    print(f"The modified DNA sequence that codes for amino acids is: 5' {target_sequence} 3'")
    
    oligo_sequence = reverse_complement(target_sequence)
    print(f"The oligo must be the reverse complement of this target.")
    print(f"\nFinal oligo sequence: 5' {oligo_sequence} 3'")
    
    return oligo_sequence

# Run the solver and print the final answer in the required format
final_answer = solve_oligo_puzzle()
print(f"\n<<<{final_answer}>>>")
