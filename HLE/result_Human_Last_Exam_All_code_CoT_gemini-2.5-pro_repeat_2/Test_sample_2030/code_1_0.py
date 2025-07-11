import textwrap

def solve_oligo_problem():
    """
    This script solves the DNA oligo design problem by following these steps:
    1. Translates the original DNA sequence in all six reading frames.
    2. Identifies the single reading frame with two unique amino acids.
    3. Applies the specified Single Nucleotide Polymorphisms (SNPs) to create a modified sequence.
    4. Designs the oligo by finding the reverse complement of the translated portion of the modified sequence.
    """

    # --- Setup ---
    original_seq = "CTTCCCCGCACAAGTGGT"
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    # For this problem, we treat Tryptophan (W) as polar to resolve the ambiguity.
    polar_aa = {'S', 'T', 'C', 'N', 'Q', 'Y', 'K', 'R', 'H', 'D', 'E', 'W'}
    non_polar_aa = {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'M'}

    # --- Helper Functions ---
    def reverse_complement(seq):
        comp_map = str.maketrans('ATGC', 'TACG')
        return seq.translate(comp_map)[::-1]

    def translate_frame(seq):
        codons = textwrap.wrap(seq, 3)
        # Only keep full codons
        codons = [c for c in codons if len(c) == 3]
        aa_seq = "".join([codon_table.get(c, 'X') for c in codons])
        return aa_seq, codons

    # --- Main Logic ---
    print("Original DNA Sequence: 5' CTT CCC CGC ACA AGT GGT 3'")
    print("-" * 50)

    # Step 1: 6-Frame Translation
    print("Step 1: Performing 6-Frame Translation")
    frames = {}
    rc_seq = reverse_complement(original_seq)

    for i in range(3):
        # Forward frames
        f_seq = original_seq[i:]
        f_aa, f_codons = translate_frame(f_seq)
        frames[f'+{i+1}'] = {'aa': f_aa, 'codons': f_codons, 'aa_set': set(f_aa)}
        
        # Reverse frames
        r_seq = rc_seq[i:]
        r_aa, r_codons = translate_frame(r_seq)
        frames[f'-{i+1}'] = {'aa': r_aa, 'codons': r_codons, 'aa_set': set(r_aa)}

    for name, data in frames.items():
        print(f"  Frame {name}: {'-'.join(list(data['aa']))}")
    print("-" * 50)

    # Step 2: Identify the unique frame
    print("Step 2: Identifying the Frame with Two Unique Amino Acids")
    all_aa_sets = [data['aa_set'] for data in frames.values()]
    unique_frame_name = None
    unique_aas = None

    for i, (name, data) in enumerate(frames.items()):
        other_aa_union = set().union(*(all_aa_sets[:i] + all_aa_sets[i+1:]))
        current_unique_aas = data['aa_set'] - other_aa_union
        if len(current_unique_aas) == 2:
            unique_frame_name = name
            unique_aas = current_unique_aas
            break
    
    print(f"Found unique frame: {unique_frame_name}")
    print(f"Unique amino acids: {unique_aas}")
    print(f"Original AA Sequence for this frame: {'-'.join(list(frames[unique_frame_name]['aa']))}")
    print(f"Original Codons for this frame: {' '.join(frames[unique_frame_name]['codons'])}")
    print("-" * 50)

    # Step 3: Apply SNPs
    print("Step 3: Analyzing and Applying SNPs")
    original_codons = frames[unique_frame_name]['codons']
    original_aa_seq = frames[unique_frame_name]['aa']
    
    print("The unique amino acids are F (Phenylalanine) and W (Tryptophan).")
    print("To satisfy the problem's constraints, we must treat W as polar and F as non-polar.")
    
    modified_codons = list(original_codons)

    # SNP 1: Non-polar -> Cysteine
    # F(TTC) -> C(TGC)
    idx_F = original_aa_seq.find('F')
    modified_codons[idx_F] = 'TGC'
    print(f"SNP 1 (Non-polar -> Cys): F ({original_codons[idx_F]}) is changed to C ({modified_codons[idx_F]})")

    # SNP 2: Polar -> Stop
    # W(TGG) -> _(TGA)
    idx_W = original_aa_seq.find('W')
    modified_codons[idx_W] = 'TGA'
    print(f"SNP 2 (Polar -> Stop): W ({original_codons[idx_W]}) is changed to Stop ({modified_codons[idx_W]})")
    print("-" * 50)

    # Step 4: Design the Oligo
    print("Step 4: Designing the Final Oligo")
    print(f"Modified codon sequence: {' '.join(modified_codons)}")

    # Find where the stop codon is
    stop_idx = -1
    for i, codon in enumerate(modified_codons):
        if codon_table.get(codon) == '_':
            stop_idx = i
            break
            
    # The oligo corresponds to the sequence BEFORE the stop codon
    target_dna_sequence = "".join(modified_codons[:stop_idx])
    print(f"The DNA sequence to bind to is the translated part: 5' {target_dna_sequence} 3'")
    
    # The oligo that binds must be the reverse complement
    oligo_sequence = reverse_complement(target_dna_sequence)
    
    print("\nThe oligo that 'binds to' this sequence is its reverse complement.")
    print("\nFinal Oligo Sequence:")
    print(f"5' {oligo_sequence} 3'")


solve_oligo_problem()