import collections

def solve_oligo_puzzle():
    """
    This function solves the bioinformatics puzzle by following the described plan.
    """
    # Step 1: Define data structures
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

    amino_acid_properties = {
        'S': 'polar', 'T': 'polar', 'C': 'polar', 'N': 'polar', 'Q': 'polar', 'Y': 'polar',
        'R': 'polar', 'H': 'polar', 'K': 'polar', 'D': 'polar', 'E': 'polar',
        'G': 'non-polar', 'A': 'non-polar', 'V': 'non-polar', 'L': 'non-polar',
        'I': 'non-polar', 'P': 'non-polar', 'F': 'non-polar', 'M': 'non-polar', 'W': 'non-polar',
    }

    # Helper functions
    def get_reverse_complement(dna):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement.get(base, base) for base in reversed(dna))

    def translate_frame(dna_frame):
        codons = [dna_frame[i:i+3] for i in range(0, len(dna_frame) - len(dna_frame) % 3, 3)]
        aa_seq = "".join([codon_table.get(c, 'X') for c in codons])
        return aa_seq, codons

    # Step 1 & 2: Translate all frames and identify the unique one
    frames = []
    # Forward frames
    for i in range(3):
        frames.append({'id': f'F{i+1}', 'dna': original_seq[i:]})
    # Reverse frames
    rev_comp = get_reverse_complement(original_seq)
    for i in range(3):
        frames.append({'id': f'F{i+4}', 'dna': rev_comp[i:]})

    all_aas_by_frame = {}
    for frame in frames:
        aa_seq, codon_seq = translate_frame(frame['dna'])
        frame['aa_seq'] = aa_seq
        frame['codons'] = codon_seq
        all_aas_by_frame[frame['id']] = set(aa_seq)

    target_frame_data = None
    for frame_id, frame_aas in all_aas_by_frame.items():
        other_aas = set()
        for other_id, other_frame_aas in all_aas_by_frame.items():
            if frame_id != other_id:
                other_aas.update(other_frame_aas)
        
        unique_aas = frame_aas - other_aas
        if len(unique_aas) == 2:
            target_frame_data = next(f for f in frames if f['id'] == frame_id)
            break
            
    # Step 3: Characterize unique amino acids
    unique_aa_info = []
    other_frames_aas = set().union(*(v for k, v in all_aas_by_frame.items() if k != target_frame_data['id']))
    for i, aa in enumerate(target_frame_data['aa_seq']):
        if aa not in other_frames_aas:
            unique_aa_info.append({
                'aa': aa,
                'codon': target_frame_data['codons'][i],
                'property': amino_acid_properties[aa],
            })

    # Step 4: Simulate SNPs, preferring transitions
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    modified_codons = {}

    for info in unique_aa_info:
        original_codon = info['codon']
        found_snps = {'transition': None, 'transversion': None}

        for i in range(3):
            original_base = original_codon[i]
            for new_base in "ACGT":
                if original_base == new_base: continue
                
                new_codon_list = list(original_codon)
                new_codon_list[i] = new_base
                new_codon = "".join(new_codon_list)
                new_aa = codon_table.get(new_codon, 'X')

                # Condition 1: Polar to Stop
                if info['property'] == 'polar' and new_aa == '_':
                    is_transition = (original_base in purines and new_base in purines) or \
                                    (original_base in pyrimidines and new_base in pyrimidines)
                    found_snps['transition' if is_transition else 'transversion'] = new_codon

                # Condition 2: Non-polar to Cysteine
                if info['property'] == 'non-polar' and new_aa == 'C':
                    is_transition = (original_base in purines and new_base in purines) or \
                                    (original_base in pyrimidines and new_base in pyrimidines)
                    found_snps['transition' if is_transition else 'transversion'] = new_codon
        
        modified_codons[original_codon] = found_snps['transition'] or found_snps['transversion']

    # Step 5: Determine the modified sequence up to the stop codon
    final_codons = list(target_frame_data['codons'])
    for i, codon in enumerate(final_codons):
        if codon in modified_codons:
            final_codons[i] = modified_codons[codon]

    translated_codons = []
    for codon in final_codons:
        if codon_table[codon] == '_':
            break
        translated_codons.append(codon)
    
    modified_translated_dna = "".join(translated_codons)

    # Step 6: Design the oligo
    oligo_sequence = get_reverse_complement(modified_translated_dna)

    # Print the final result
    print(f"The DNA sequence for the oligo is: 5' {oligo_sequence} 3'")

solve_oligo_puzzle()