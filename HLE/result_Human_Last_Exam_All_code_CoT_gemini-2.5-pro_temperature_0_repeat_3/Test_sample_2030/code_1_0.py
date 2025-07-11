import textwrap

def solve_oligo_design():
    """
    Solves the DNA oligo design problem based on translation, SNPs, and binding.
    """
    # --- Step 1: Define constants and helper functions ---
    original_seq = "CTTCCCGCACAAGTGGT"
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
    # To solve the problem, we must assume Tryptophan (W) is the polar AA and Phenylalanine (F) is the non-polar AA.
    polar_aa = {'S', 'T', 'N', 'Q', 'Y', 'K', 'R', 'H', 'D', 'E', 'W'}
    non_polar_aa = {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'M'}

    def reverse_complement(seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement.get(base, base) for base in reversed(seq))

    def translate_frame(seq):
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        amino_acids = []
        dna_codons = []
        for codon in codons:
            if len(codon) == 3:
                amino_acids.append(codon_table.get(codon, '?'))
                dna_codons.append(codon)
        return amino_acids, dna_codons

    print(f"Original Sequence: 5' {original_seq} 3'")
    
    # --- Step 2: Generate and translate all 6 frames ---
    translations = []
    frame_codons = []
    # Forward frames
    for i in range(3):
        aa_seq, codon_seq = translate_frame(original_seq[i:])
        translations.append(aa_seq)
        frame_codons.append(codon_seq)
    # Reverse complement frames
    rev_comp_seq = reverse_complement(original_seq)
    for i in range(3):
        aa_seq, codon_seq = translate_frame(rev_comp_seq[i:])
        translations.append(aa_seq)
        frame_codons.append(codon_seq)

    # --- Step 3: Identify the unique frame ---
    aa_counts_per_frame = [set(t) for t in translations]
    unique_frame_index = -1
    unique_aas_in_frame = []
    for i, current_frame_aas in enumerate(aa_counts_per_frame):
        other_frames_aas = set()
        for j, other_frame_aas in enumerate(aa_counts_per_frame):
            if i != j:
                other_frames_aas.update(other_frame_aas)
        unique_to_this_frame = current_frame_aas - other_frames_aas
        if len(unique_to_this_frame) == 2:
            unique_frame_index = i
            unique_aas_in_frame = list(unique_to_this_frame)
            break
    
    print(f"\nIdentified reading frame number {unique_frame_index + 1} as the correct one.")
    print(f"It contains two unique amino acids: {', '.join(unique_aas_in_frame)}.")

    # --- Step 4: Identify target amino acids and their codons ---
    target_frame_translation = translations[unique_frame_index]
    target_frame_codons = frame_codons[unique_frame_index]
    polar_target_codon, non_polar_target_codon = None, None
    polar_target_aa, non_polar_target_aa = None, None

    for aa, codon in zip(target_frame_translation, target_frame_codons):
        if aa in unique_aas_in_frame:
            if aa in polar_aa:
                polar_target_aa, polar_target_codon = aa, codon
            elif aa in non_polar_aa:
                non_polar_target_aa, non_polar_target_codon = aa, codon
    
    print(f"The polar amino acid is {polar_target_aa} (from codon {polar_target_codon}).")
    print(f"The non-polar amino acid is {non_polar_target_aa} (from codon {non_polar_target_codon}).")

    # --- Step 5: Simulate the SNPs ---
    # Non-polar to Cysteine
    cysteine_codons = {'TGT', 'TGC'}
    modified_non_polar_codon = next(c for c in cysteine_codons if sum(1 for a, b in zip(non_polar_target_codon, c) if a != b) == 1)
    # Polar to Stop
    stop_codons = {'TAA', 'TAG', 'TGA'}
    modified_polar_codon = next(c for c in stop_codons if sum(1 for a, b in zip(polar_target_codon, c) if a != b) == 1)

    print(f"\nSNP 1: Non-polar {non_polar_target_aa} ({non_polar_target_codon}) -> Cysteine ({modified_non_polar_codon})")
    print(f"SNP 2: Polar {polar_target_aa} ({polar_target_codon}) -> Stop ({modified_polar_codon})")

    # --- Step 6: Construct the modified sequence and design the oligo ---
    modified_codons = [
        modified_non_polar_codon if c == non_polar_target_codon else
        modified_polar_codon if c == polar_target_codon else c
        for c in target_frame_codons
    ]
    
    oligo_template_codons = []
    for codon in modified_codons:
        if codon_table[codon] == '_':
            break
        oligo_template_codons.append(codon)

    oligo_template_sequence = "".join(oligo_template_codons)
    oligo_sequence = reverse_complement(oligo_template_sequence)

    print(f"\nThe modified DNA sequence that is translated into amino acids is: 5' {' '.join(textwrap.wrap(oligo_template_sequence, 3))} 3'")
    print("\nThe oligo must be the reverse complement of this sequence to bind to it.")
    print("\nFinal Oligo Sequence:")
    # The final "equation" is the sequence itself. We print each base.
    final_output_string = f"5' {' '.join(list(oligo_sequence))} 3'"
    print(final_output_string)
    
    # Final answer in specified format
    print(f"\n<<<{oligo_sequence}>>>")

solve_oligo_design()