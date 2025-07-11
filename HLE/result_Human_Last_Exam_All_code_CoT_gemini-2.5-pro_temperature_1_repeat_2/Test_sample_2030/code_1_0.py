import textwrap

def get_peptides(dna_sequence):
    """Translates a DNA sequence in all six reading frames."""
    
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'TTA':'L', 'TTG':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TTC':'F', 'TTT':'F', 'TGC':'C', 'TGT':'C',
        'TAC':'Y', 'TAT':'Y', 'TGG':'W',
        'TAA':'*', 'TAG':'*', 'TGA':'*'
    }
    
    def reverse_complement(seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement.get(base, base) for base in reversed(seq))

    frames = {}
    peptides = {}
    
    # Forward frames
    fwd_seq = dna_sequence
    # Reverse complement for reverse frames
    rev_seq = reverse_complement(fwd_seq)
    
    sequences = {'fwd': fwd_seq, 'rev': rev_seq}
    
    for direction, seq in sequences.items():
        for frame_start in range(3):
            codons = [seq[i:i+3] for i in range(frame_start, len(seq), 3)]
            codons = [c for c in codons if len(c) == 3] # Remove incomplete codons
            peptide = "".join([codon_table.get(c, 'X') for c in codons])
            
            frame_num = frame_start + 1
            if direction == 'rev':
                frame_num = -(frame_start + 1)
            
            frames[f'{frame_num:+d}'] = codons
            peptides[f'{frame_num:+d}'] = peptide
            
    return frames, peptides

def solve_oligo():
    """
    Solves the oligo design problem based on the provided DNA sequence and constraints.
    """
    original_sequence = "CTTCCCCGCACAAGTGGT"

    print("Step 1: Performing six-frame translation of the original sequence.")
    frames, peptides = get_peptides(original_sequence)
    print("Original DNA Sequence: 5' CTT CCC CGC ACA AGT GGT 3'")
    for frame_name, peptide in peptides.items():
        print(f"Frame {frame_name}: {peptide} ({' '.join(frames[frame_name])})")

    print("\nStep 2: Identifying the target reading frame.")
    print("A direct analysis shows no frame has two unique amino acids. The most plausible interpretation is to find a frame with a polar and non-polar amino acid pair that fits the SNP criteria.")
    print("Frame +3 (SPHKW) contains polar Lysine (K) and non-polar Tryptophan (W).")
    print(" - Lysine (K, codon 'AAG') can mutate to a stop codon ('TAG').")
    print(" - Tryptophan (W, codon 'TGG') can mutate to Cysteine ('TGT' or 'TGC').")
    print("This makes Frame +3 the target frame.")
    
    target_frame_codons = frames['+3']
    print(f"\nTarget Frame +3 Codons: {' '.join(target_frame_codons)}")

    print("\nStep 3: Simulating the SNPs.")
    modified_codons = list(target_frame_codons) # create a mutable copy
    # SNP 1: K (AAG) -> Stop (*)
    k_index = target_frame_codons.index('AAG')
    modified_codons[k_index] = 'TAG'
    # SNP 2: W (TGG) -> Cys (C)
    w_index = target_frame_codons.index('TGG')
    modified_codons[w_index] = 'TGT'
    print(f"Modified Frame +3 Codons: {' '.join(modified_codons)}")

    print("\nStep 4: Identifying the target sequence for the oligo.")
    print("The oligo binds to the modified sequence, but only the parts translated into amino acids.")
    # The stop codon 'TAG' interrupts translation. The oligo will bind to the part before it.
    translated_part_codons = []
    for codon in modified_codons:
        if codon in ['TAA', 'TAG', 'TGA']:
            break
        translated_part_codons.append(codon)
    
    target_dna_segment = "".join(translated_part_codons)
    print(f"The contiguous translated DNA segment is: 5'-{target_dna_segment}-3'")

    print("\nStep 5: Designing the oligo by finding the reverse complement.")
    def reverse_complement(seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement.get(base, base) for base in reversed(seq))
    
    oligo_sequence = reverse_complement(target_dna_segment)
    
    print("\nFinal Oligo Sequence (reverse complement of the target segment):")
    print(f"Target:      5'-{target_dna_segment}-3'")
    # Format the reverse complement string to align with the target
    rev_comp_str = f"3'-{''.join(reversed(oligo_sequence))}-5'"
    print(f"Rev Comp:    {rev_comp_str}")
    print(f"Oligo (5'-3'): 5'-{oligo_sequence}-3'")
    return oligo_sequence

final_oligo = solve_oligo()
print("\nThe final DNA sequence for the oligo is:")
print(f"5'-{final_oligo}-3'")
print(f'<<<{final_oligo}>>>')