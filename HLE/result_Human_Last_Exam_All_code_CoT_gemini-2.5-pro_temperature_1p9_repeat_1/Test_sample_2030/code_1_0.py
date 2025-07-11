def solve_oligo_puzzle():
    """
    Solves the DNA puzzle to find the required oligo sequence.
    """
    original_seq = "CTTCCCCGCACAAGTGGT"

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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    # Polar: S, T, N, Q, Y. Non-polar: A, V, L, I, P, F, W, M, G. (C is target, so we classify the original AA)
    amino_acid_properties = {
        'A': 'non-polar', 'V': 'non-polar', 'L': 'non-polar', 'I': 'non-polar',
        'P': 'non-polar', 'F': 'non-polar', 'W': 'non-polar', 'M': 'non-polar', 'G': 'non-polar',
        'S': 'polar', 'T': 'polar', 'N': 'polar', 'Q': 'polar', 'Y': 'polar', 'C': 'polar'
    }

    def get_reverse_complement(dna_seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement[base] for base in reversed(dna_seq))

    def translate_frame(dna_seq):
        codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
        aa_seq = ""
        codon_list = []
        for codon in codons:
            if len(codon) == 3:
                aa_seq += codon_table.get(codon, '?')
                codon_list.append(codon)
        return aa_seq, codon_list

    print("Step 1: Translating all six reading frames...")
    
    # Forward frames
    frame1_aa, frame1_codons = translate_frame(original_seq)
    frame2_aa, frame2_codons = translate_frame(original_seq[1:])
    frame3_aa, frame3_codons = translate_frame(original_seq[2:])
    
    # Reverse frames
    rev_comp_seq = get_reverse_complement(original_seq)
    frame4_aa, frame4_codons = translate_frame(rev_comp_seq)
    frame5_aa, frame5_codons = translate_frame(rev_comp_seq[1:])
    frame6_aa, frame6_codons = translate_frome(rev_comp_seq[2:])

    all_frames = {
        "+1": (frame1_aa, frame1_codons), "+2": (frame2_aa, frame2_codons), "+3": (frame3_aa, frame3_codons),
        "-1": (frame4_aa, frame4_codons), "-2": (frame5_aa, frame5_codons), "-3": (frame6_aa, frame6_codons)
    }
    for name, (aa_seq, _) in all_frames.items():
        print(f"Frame {name}: {aa_seq}")
        
    print("\nStep 2: Finding the reading frame with two unique amino acids...")
    
    all_aa = "".join([f[0] for f in all_frames.values()])
    unique_frame_name = None
    unique_aa_info = []

    for name, (aa_seq, codons) in all_frames.items():
        other_aa = "".join([f[0] for fname, f in all_frames.items() if fname != name])
        current_unique_aa = []
        for i, aa in enumerate(aa_seq):
            if aa in amino_acid_properties and other_aa.count(aa) == 0:
                 current_unique_aa.append({'aa': aa, 'codon': codons[i], 'prop': amino_acid_properties[aa]})
        
        if len(current_unique_aa) == 2:
            unique_frame_name = name
            unique_aa_info = current_unique_aa
            print(f"Found unique frame: {name}")
            print(f"Amino Acid Sequence: {aa_seq}")
            print(f"The two unique amino acids are {unique_aa_info[0]['aa']} ({unique_aa_info[0]['codon']}) and {unique_aa_info[1]['aa']} ({unique_aa_info[1]['codon']}).")
            break
            
    print("\nStep 3: Applying SNPs to the unique codons...")
    
    # The problem logic: polar -> stop, non-polar -> Cys
    modified_codons = {}
    for info in unique_aa_info:
        if info['prop'] == 'polar':
            # Find SNP for polar -> stop
            original_codon = info['codon']
            for i in range(3):
                for base in "ATGC":
                    if base != original_codon[i]:
                        new_codon = original_codon[:i] + base + original_codon[i+1:]
                        if codon_table.get(new_codon) == '_':
                            modified_codons['stop'] = new_codon
                            print(f"SNP found: {info['aa']} ({original_codon}) -> Stop ({new_codon})")
                            break
                if 'stop' in modified_codons: break
        
        elif info['prop'] == 'non-polar':
            # Find SNP for non-polar -> Cys
            original_codon = info['codon']
            for i in range(3):
                for base in "ATGC":
                    if base != original_codon[i]:
                        new_codon = original_codon[:i] + base + original_codon[i+1:]
                        if codon_table.get(new_codon) == 'C':
                            modified_codons['cys'] = new_codon
                            print(f"SNP found: {info['aa']} ({original_codon}) -> Cysteine ({new_codon})")
                            break
                if 'cys' in modified_codons: break
    
    print("\nStep 4: Constructing the oligo target sequence...")

    original_codons_in_frame = all_frames[unique_frame_name][1]
    
    # Identify which original codons were unique to build the target
    original_unique_codons = [info['codon'] for info in unique_aa_info]

    # Reconstruct the sequence up to the stop codon
    target_dna_sequence = ""
    stop_codon_created = False
    for codon in original_codons_in_frame:
        if codon in original_unique_codons:
            if codon == [info['codon'] for info in unique_aa_info if info['prop'] == 'polar'][0]:
                print(f"The oligo sequence will not include the sequence for the new stop codon ({modified_codons['stop']}) or any downstream sequences.")
                stop_codon_created = True
                break
            else: # It must be the non-polar one
                target_dna_sequence += modified_codons['cys']
        else:
            target_dna_sequence += codon

    print(f"The DNA sequence translated into amino acids in the modified frame is: 5' {target_dna_sequence} 3'")
    
    print("\nStep 5: Calculating the final oligo sequence (reverse complement of the target)...")
    
    oligo_sequence = get_reverse_complement(target_dna_sequence)
    
    print("\nThe final oligo is the reverse complement of the target DNA.")
    print(f"Final Oligo Sequence: 5' {oligo_sequence} 3'")
    
    return oligo_sequence

final_answer = solve_oligo_puzzle()
print(f"<<<{final_answer}>>>")