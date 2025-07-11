import textwrap

def get_peptide_and_codons(dna_sequence):
    """Translates a DNA sequence into a peptide sequence."""
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    peptides = []
    codons = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        codons.append(codon)
        peptides.append(codon_table.get(codon, 'X'))
    return "".join(peptides), codons

def reverse_complement(dna_sequence):
    """Calculates the reverse complement of a DNA sequence."""
    complement_map = str.maketrans('ATCG', 'TAGC')
    return dna_sequence.upper().translate(complement_map)[::-1]

def find_snp_codon(original_codon, target_codons):
    """Finds a single nucleotide change to convert a codon to one of the target codons."""
    bases = "ATCG"
    for i in range(3):
        for base in bases:
            if base != original_codon[i]:
                new_codon = list(original_codon)
                new_codon[i] = base
                new_codon = "".join(new_codon)
                if new_codon in target_codons:
                    return new_codon
    return None

def solve_oligo_puzzle():
    """Main function to solve the puzzle."""
    original_seq = "CTTCCCCGCACAAGTGGT"
    aa_properties = {
        'G': 'Non-polar', 'A': 'Non-polar', 'V': 'Non-polar', 'L': 'Non-polar',
        'I': 'Non-polar', 'P': 'Non-polar', 'F': 'Non-polar', 'M': 'Non-polar', 'W': 'Non-polar',
        'S': 'Polar', 'T': 'Polar', 'C': 'Polar', 'N': 'Polar', 'Q': 'Polar', 'Y': 'Polar',
        'D': 'Polar', 'E': 'Polar', 'K': 'Polar', 'R': 'Polar', 'H': 'Polar'
    }
    stop_codons = ["TAA", "TAG", "TGA"]
    cys_codons = ["TGT", "TGC"]

    print("--- Step 1: Translating all 6 reading frames ---")
    rev_comp_seq = reverse_complement(original_seq)
    
    frames_data = {}
    # Forward frames
    for frame in range(3):
        seq_to_translate = original_seq[frame:]
        peptides, codons = get_peptide_and_codons(seq_to_translate)
        frames_data[f"Forward Frame {frame+1}"] = {'peptides': peptides, 'codons': codons}
    
    # Reverse frames
    for frame in range(3):
        seq_to_translate = rev_comp_seq[frame:]
        peptides, codons = get_peptide_and_codons(seq_to_translate)
        frames_data[f"Reverse Frame {frame+1}"] = {'peptides': peptides, 'codons': codons}
        
    for name, data in frames_data.items():
        print(f"{name}: {data['peptides']} ({' '.join(data['codons'])})")
    
    print("\n--- Step 2: Identifying the frame with two unique amino acids ---")
    
    all_peptides_pool = set("".join(d['peptides'] for d in frames_data.values()))
    target_frame_name = None
    target_frame_data = None
    unique_aas = []

    for name, data in frames_data.items():
        other_peptides = set("".join(d['peptides'] for n, d in frames_data.items() if n != name))
        current_unique = [p for p in set(data['peptides']) if p not in other_peptides]
        if len(current_unique) == 2:
            target_frame_name = name
            target_frame_data = data
            unique_aas = current_unique
            break
            
    print(f"Found target frame: {target_frame_name}")
    print(f"Unique amino acids in this frame: {unique_aas[0]} and {unique_aas[1]}")

    print("\n--- Step 3: Analyzing the SNPs ---")
    
    aa1, aa2 = unique_aas[0], unique_aas[1]
    
    # Determine which is polar and which is non-polar
    if aa_properties[aa1] == 'Polar':
        polar_aa, nonpolar_aa = aa1, aa2
    else:
        polar_aa, nonpolar_aa = aa2, aa1
        
    # Get their original codons
    polar_codon = target_frame_data['codons'][target_frame_data['peptides'].find(polar_aa)]
    nonpolar_codon = target_frame_data['codons'][target_frame_data['peptides'].find(nonpolar_aa)]
    
    print(f"Polar amino acid: {polar_aa} (from codon {polar_codon})")
    print(f"Non-polar amino acid: {nonpolar_aa} (from codon {nonpolar_codon})")
    
    # Find the modified codons
    new_stop_codon = find_snp_codon(polar_codon, stop_codons)
    new_cys_codon = find_snp_codon(nonpolar_codon, cys_codons)

    print(f"SNP 1: {polar_aa} ({polar_codon}) -> Stop ({new_stop_codon})")
    print(f"SNP 2: {nonpolar_aa} ({nonpolar_codon}) -> Cys ({new_cys_codon})")

    print("\n--- Step 4: Constructing the oligo target sequence ---")
    
    modified_codons_for_oligo = []
    final_peptide_seq = []
    
    for peptide, codon in zip(target_frame_data['peptides'], target_frame_data['codons']):
        if peptide == polar_aa:
            final_peptide_seq.append("Stop")
            break # Stop translation
        elif peptide == nonpolar_aa:
            modified_codons_for_oligo.append(new_cys_codon)
            final_peptide_seq.append("Cys")
        else:
            modified_codons_for_oligo.append(codon)
            final_peptide_seq.append(peptide)

    oligo_target_dna = "".join(modified_codons_for_oligo)
    
    print(f"Modified peptide sequence: {'-'.join(final_peptide_seq)}")
    print(f"DNA sequence on template strand coding for these amino acids: 5'-{oligo_target_dna}-3'")
    
    print("\n--- Step 5: Designing the final oligo ---")
    
    final_oligo_seq = reverse_complement(oligo_target_dna)
    print("The oligo must be the reverse complement of the target DNA sequence.")
    print(f"Final Oligo Sequence: 5'-{final_oligo_seq}-3'")
    
    return final_oligo_seq

# Run the solver and print the final answer in the required format
final_answer = solve_oligo_puzzle()
print(f"\n<<<5'-{final_answer}-3'>>>")
