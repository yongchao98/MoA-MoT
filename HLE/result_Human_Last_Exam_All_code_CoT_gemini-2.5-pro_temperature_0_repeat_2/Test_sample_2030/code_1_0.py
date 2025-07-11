import textwrap

def solve_oligo_design():
    """
    Solves the DNA oligo design problem based on sequence translation and SNP analysis.
    """
    # Step 1: Define constants and helper functions
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

    polar_aa = {'S', 'T', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H'}
    non_polar_aa = {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M'}
    cysteine_codons = {'TGT', 'TGC'}
    stop_codons = {'TAA', 'TAG', 'TGA'}

    def get_reverse_complement(dna):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement.get(base, 'N') for base in reversed(dna))

    def translate_frame(dna):
        codons = textwrap.wrap(dna, 3)
        peptide = []
        for codon in codons:
            if len(codon) == 3:
                peptide.append(codon_table.get(codon, 'X'))
        return peptide

    # Step 2: Generate and translate all 6 frames
    frames_dna = []
    frames_dna.append(original_seq[0:])
    frames_dna.append(original_seq[1:])
    frames_dna.append(original_seq[2:])
    
    rev_comp = get_reverse_complement(original_seq)
    frames_dna.append(rev_comp[0:])
    frames_dna.append(rev_comp[1:])
    frames_dna.append(rev_comp[2:])

    frames_peptide = [translate_frame(dna) for dna in frames_dna]
    
    # Step 3: Identify the unique reading frame
    all_peptides_set = set()
    for i, p_list in enumerate(frames_peptide):
        all_peptides_set.update(p_list)

    unique_frame_idx = -1
    for i, p_list in enumerate(frames_peptide):
        other_peptides = set()
        for j, other_p_list in enumerate(frames_peptide):
            if i != j:
                other_peptides.update(other_p_list)
        
        current_frame_unique_aas = set(p_list) - other_peptides
        if len(current_frame_unique_aas) == 2:
            # Check if this frame has both polar and non-polar AAs for the SNPs
            has_polar = any(aa in polar_aa for aa in p_list)
            if has_polar:
                unique_frame_idx = i
                break
    
    print(f"Identified Reading Frame: Frame {unique_frame_idx + 1}")
    
    target_frame_dna_str = frames_dna[unique_frame_idx]
    target_frame_codons = textwrap.wrap(target_frame_dna_str, 3)
    target_frame_codons = [c for c in target_frame_codons if len(c) == 3] # remove incomplete codons
    target_frame_peptide = frames_peptide[unique_frame_idx]

    print(f"Original DNA in Frame: {' '.join(target_frame_codons)}")
    print(f"Original Peptide: {' '.join(target_frame_peptide)}\n")

    # Step 4: Find the SNPs
    modified_codons = list(target_frame_codons)
    stop_pos = -1
    
    # SNP 1: Polar -> Stop
    for i, aa in enumerate(target_frame_peptide):
        if aa in polar_aa:
            original_codon = target_frame_codons[i]
            for pos in range(3):
                for base in 'ATCG':
                    if base != original_codon[pos]:
                        mutated_codon = original_codon[:pos] + base + original_codon[pos+1:]
                        if mutated_codon in stop_codons:
                            print(f"SNP 1 (Polar -> Stop): Found in {aa} (codon {original_codon})")
                            print(f"Modified Codon: {mutated_codon} (Stop)")
                            modified_codons[i] = mutated_codon
                            stop_pos = i
                            break
                if stop_pos != -1: break
        if stop_pos != -1: break
    
    # SNP 2: Non-polar -> Cysteine
    cys_pos = -1
    for i, aa in enumerate(target_frame_peptide):
        # The oligo targets the sequence before the stop codon, so the SNP must be in that region
        if i < stop_pos and aa in non_polar_aa:
            original_codon = target_frame_codons[i]
            for pos in range(3):
                for base in 'ATCG':
                    if base != original_codon[pos]:
                        mutated_codon = original_codon[:pos] + base + original_codon[pos+1:]
                        if mutated_codon in cysteine_codons:
                            print(f"SNP 2 (Non-polar -> Cys): Found in {aa} (codon {original_codon})")
                            print(f"Modified Codon: {mutated_codon} (Cys)\n")
                            modified_codons[i] = mutated_codon
                            cys_pos = i
                            break
                if cys_pos != -1: break
        if cys_pos != -1: break

    # Step 5: Design the oligo
    print(f"Modified DNA in Frame: {' '.join(modified_codons)}")
    
    # The oligo binds to the part of the sequence translated into amino acids
    oligo_target_dna = "".join(modified_codons[:stop_pos])
    print(f"Oligo Target Sequence (pre-stop codon): 5' {oligo_target_dna} 3'")
    
    final_oligo_seq = get_reverse_complement(oligo_target_dna)
    print(f"\nFinal Oligo Sequence (Reverse Complement): 5' {final_oligo_seq} 3'")
    
    return final_oligo_seq

# Execute the function and print the final answer in the required format.
final_answer = solve_oligo_design()
print(f"\n<<<{final_answer}>>>")
