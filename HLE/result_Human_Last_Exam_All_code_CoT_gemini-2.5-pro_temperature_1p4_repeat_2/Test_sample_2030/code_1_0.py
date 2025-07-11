import collections

def solve_oligo_puzzle():
    """
    Solves the DNA oligo design puzzle by translating, identifying a unique frame,
    applying specified SNPs, and designing a complementary oligo.
    """

    # Step 1: Establish Genetic Data
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    stop_codons = {'TAA', 'TAG', 'TGA'}
    cys_codons = {'TGC', 'TGT'}

    polar_aa = {'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H'}
    non_polar_aa = {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W'}

    original_seq = "CTTCCCGCGACAAGTGGT"

    def get_reverse_complement(dna):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement.get(base, base) for base in reversed(dna))

    # Step 2: Generate Six-Frame Translation
    frames = []
    # Forward frames
    for i in range(3):
        codons = [original_seq[j:j+3] for j in range(i, len(original_seq), 3) if j+3 <= len(original_seq)]
        peptide = "".join([codon_table.get(c, 'X') for c in codons])
        frames.append({'frame': i + 1, 'codons': codons, 'peptide': peptide})
    
    # Reverse frames
    rev_comp_seq = get_reverse_complement(original_seq)
    for i in range(3):
        codons = [rev_comp_seq[j:j+3] for j in range(i, len(rev_comp_seq), 3) if j+3 <= len(rev_comp_seq)]
        peptide = "".join([codon_table.get(c, 'X') for c in codons])
        frames.append({'frame': -(i + 1), 'codons': codons, 'peptide': peptide})

    print("--- Identifying Unique Reading Frame ---")
    all_peptides = [set(f['peptide']) for f in frames]
    unique_frame_data = None
    for i in range(len(frames)):
        other_aas = set()
        for j in range(len(frames)):
            if i == j:
                continue
            other_aas.update(all_peptides[j])
        
        current_aas = all_peptides[i]
        unique_aas = current_aas - other_aas
        
        if len(unique_aas) == 2:
            unique_frame_data = frames[i]
            print(f"Found unique frame: Frame {unique_frame_data['frame']}")
            print(f"Sequence: {' '.join(unique_frame_data['codons'])}")
            print(f"Peptide: {unique_frame_data['peptide']}")
            print(f"Unique amino acids: {', '.join(unique_aas)}\n")
            break

    # Step 3 & 4: Pinpoint SNP Locations
    print("--- Analyzing SNPs for the Modified Sequence ---")
    codons_for_aa = collections.defaultdict(list)
    for codon, aa in codon_table.items():
        codons_for_aa[aa].append(codon)

    def is_single_snp(codon1, codon2):
        diffs = 0
        for c1, c2 in zip(codon1, codon2):
            if c1 != c2:
                diffs += 1
        return diffs == 1

    original_codons = unique_frame_data['codons']
    original_peptide = unique_frame_data['peptide']
    modified_codons = list(original_codons)

    snp1_found = False
    snp2_found = False

    for i, (codon, aa) in enumerate(zip(original_codons, original_peptide)):
        # Condition 1: Polar amino acid to a stop codon
        if aa in polar_aa:
            for stop in stop_codons:
                if is_single_snp(codon, stop):
                    print(f"SNP 1: Polar amino acid '{aa}' ({codon}) changed to Stop Codon ({stop}).")
                    modified_codons[i] = stop
                    snp1_found = True
                    break
        if snp1_found: break

    for i, (codon, aa) in enumerate(zip(original_codons, original_peptide)):
        # Condition 2: Non-polar amino acid to a Cysteine
        if aa in non_polar_aa:
            for cys in cys_codons:
                if is_single_snp(codon, cys):
                    print(f"SNP 2: Non-polar amino acid '{aa}' ({codon}) changed to Cysteine ({cys}).\n")
                    modified_codons[i] = cys
                    snp2_found = True
                    break
        if snp2_found: break
    
    # Step 5 & 6: Construct the Modified Sequence and Design the Oligo
    print("--- Designing the Oligo ---")
    
    # Find the new stop codon and get the sequence before it
    translated_part = []
    for codon in modified_codons:
        if codon in stop_codons:
            break
        translated_part.append(codon)
    
    modified_target_seq = "".join(translated_part)
    print(f"Modified sequence that is translated: 5' {modified_target_seq} 3'")
    
    oligo_sequence = get_reverse_complement(modified_target_seq)
    print(f"\nFinal Oligo Sequence (reverse complement of the translated part):")
    print(f"5' {oligo_sequence} 3'")

    return oligo_sequence

final_answer = solve_oligo_puzzle()
# The final answer is wrapped for programmatic retrieval.
# The code above prints the step-by-step derivation.
print(f"\n<<<5' {final_answer} 3'>>>")
