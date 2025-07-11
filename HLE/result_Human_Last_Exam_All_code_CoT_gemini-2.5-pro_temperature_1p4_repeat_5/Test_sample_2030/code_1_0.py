def solve_oligo_design():
    """
    Solves the oligo design problem by performing a six-frame translation,
    identifying the correct reading frame based on unique amino acids,
    simulating the specified SNPs, and designing the corresponding oligo.
    """

    # --- Step 1: Define data and helper functions ---
    original_seq = "CTTCCCGCACAAGTGGT"
    
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
        'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',
        'TGT': 'C', 'TGC': 'C', 'TGA': 'Stop', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
        'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
        'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
        'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    aa_properties = {
        'A': 'nonpolar', 'V': 'nonpolar', 'I': 'nonpolar', 'L': 'nonpolar',
        'M': 'nonpolar', 'F': 'nonpolar', 'Y': 'polar', 'W': 'nonpolar',
        'S': 'polar', 'T': 'polar', 'C': 'polar', 'N': 'polar', 'Q': 'polar',
        'H': 'basic', 'K': 'basic', 'R': 'basic', 'D': 'acidic', 'E': 'acidic',
        'G': 'nonpolar', 'P': 'nonpolar'
    }

    def get_reverse_complement(dna_seq):
        complement_map = str.maketrans('ATCG', 'TAGC')
        return dna_seq.translate(complement_map)[::-1]

    # --- Step 2: Six-frame translation ---
    frames = []
    peptides = []
    
    # Forward frames
    for i in range(3):
        codons = [original_seq[j:j+3] for j in range(i, len(original_seq), 3) if len(original_seq[j:j+3]) == 3]
        frames.append({'dna_codons': codons, 'dna_seq': original_seq[i:]})
        peptides.append([codon_table.get(c, '?') for c in codons])

    # Reverse frames
    rev_comp_seq = get_reverse_complement(original_seq)
    for i in range(3):
        codons = [rev_comp_seq[j:j+3] for j in range(i, len(rev_comp_seq), 3) if len(rev_comp_seq[j:j+3]) == 3]
        frames.append({'dna_codons': codons, 'dna_seq': rev_comp_seq[i:]})
        peptides.append([codon_table.get(c, '?') for c in codons])

    # --- Step 3: Identify the target frame ---
    target_frame_idx = -1
    unique_aas_in_target = []
    
    for i in range(len(peptides)):
        current_aas = set(peptides[i])
        other_aas = set()
        for j in range(len(peptides)):
            if i != j:
                other_aas.update(peptides[j])
        
        unique_aas = current_aas - other_aas
        if len(unique_aas) == 2:
            target_frame_idx = i
            unique_aas_in_target = list(unique_aas)
            break
            
    # --- Step 4: Analyze SNPs in the target frame ---
    target_frame_codons = frames[target_frame_idx]['dna_codons']
    target_peptide = peptides[target_frame_idx]
    
    unique_codon_map = {}
    for aa, codon in zip(target_peptide, target_frame_codons):
        if aa in unique_aas_in_target:
            unique_codon_map[aa] = codon

    polar_aa = None
    nonpolar_aa = None
    for aa in unique_aas_in_target:
        if aa_properties[aa] == 'polar':
            polar_aa = aa
        elif aa_properties[aa] == 'nonpolar':
            nonpolar_aa = aa

    polar_codon = unique_codon_map[polar_aa]       # This will be 'CAA' for Q
    nonpolar_codon = unique_codon_map[nonpolar_aa] # This will be 'TTC' for F

    # One SNP changed a polar amino acid to a stop codon
    # 'CAA' (Q) -> 'TAA' (Stop) is a single C->T change (transition).
    modified_polar_codon = 'TAA'
    
    # The other SNP changed a non-polar amino acid to a cysteine
    # 'TTC' (F) -> 'TGT' (C) is a single C->T change (transition).
    # 'TTC' (F) -> 'TGC' (C) is a single T->G change (transversion).
    # We choose the transition to be consistent.
    modified_nonpolar_codon = 'TGT'
    
    # --- Step 5: Construct modified sequence and design oligo ---
    modified_codons = []
    for codon in target_frame_codons:
        if codon == polar_codon:
            modified_codons.append(modified_polar_codon)
        elif codon == nonpolar_codon:
            modified_codons.append(modified_nonpolar_codon)
        else:
            modified_codons.append(codon)
    
    # Isolate the part before the stop codon
    oligo_target_codons = []
    for codon in modified_codons:
        if codon_table.get(codon) == 'Stop':
            break
        oligo_target_codons.append(codon)
        
    oligo_target_sequence = "".join(oligo_target_codons)
    
    # The oligo sequence is the reverse complement of the target sequence
    oligo_sequence = get_reverse_complement(oligo_target_sequence)
    
    print("Original DNA Sequence: 5' CTT CCC CGC ACA AGT GGT 3'")
    print(f"Identified Target Reading Frame: {target_frame_idx + 1 if target_frame_idx < 3 else -(target_frame_idx - 2)}")
    print(f"Original Peptide: {'-'.join(target_peptide)}")
    print(f"Unique Amino Acids and their codons: {polar_aa} ({polar_codon}), {nonpolar_aa} ({nonpolar_codon})")
    print(f"SNP Changes: {polar_codon} -> {modified_polar_codon} (Stop), {nonpolar_codon} -> {modified_nonpolar_codon} (Cys)")
    print(f"Modified DNA fragment (pre-stop): 5' {' '.join(oligo_target_codons)} 3'")
    print("-" * 20)
    print(f"The DNA sequence of the oligo that binds to this modified fragment is:")
    print(f"5' {oligo_sequence} 3'")
    return oligo_sequence

final_answer = solve_oligo_design()
# The final answer is requested in the special format below.
# <<<TGCGGGCAC>>>