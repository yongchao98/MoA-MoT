import textwrap

def get_reverse_complement(dna_seq):
    """Calculates the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATGC", "TACG")
    return dna_seq.upper().translate(complement_map)[::-1]

def translate_dna(dna_seq):
    """Translates a DNA sequence into an amino acid sequence."""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    peptide = ""
    for i in range(0, len(dna_seq) - (len(dna_seq) % 3), 3):
        codon = dna_seq[i:i+3]
        peptide += codon_table.get(codon, 'X')
    return peptide

def main():
    """Main function to solve the problem."""
    # Step 1: Define initial sequences and properties
    original_seq = "CTTCCCGCACAAGTGGT"
    rev_comp_seq = get_reverse_complement(original_seq)
    
    aa_properties = {
        'F': ('Phe', 'Non-polar'), 'L': ('Leu', 'Non-polar'), 'I': ('Ile', 'Non-polar'), 
        'M': ('Met', 'Non-polar'), 'V': ('Val', 'Non-polar'), 'S': ('Ser', 'Polar'), 
        'P': ('Pro', 'Non-polar'), 'T': ('Thr', 'Polar'), 'A': ('Ala', 'Non-polar'), 
        'Y': ('Tyr', 'Polar'), 'H': ('His', 'Polar'), 'Q': ('Gln', 'Polar'), 
        'N': ('Asn', 'Polar'), 'K': ('Lys', 'Polar'), 'D': ('Asp', 'Polar'), 
        'E': ('Glu', 'Polar'), 'C': ('Cys', 'Polar'), 'W': ('Trp', 'Non-polar'), 
        'R': ('Arg', 'Polar'), 'G': ('Gly', 'Non-polar')
    }
    stop_codons = ['TAA', 'TAG', 'TGA']
    cys_codons = ['TGT', 'TGC']

    print("--- Step 1: Sequence Translation ---")
    print(f"Original Sequence (5'->3'):  {textwrap.fill(original_seq, 18)}")
    print(f"Reverse Complement (5'->3'): {textwrap.fill(rev_comp_seq, 18)}")
    print("-" * 35)

    # Step 2: Translate all 6 reading frames
    frames = {}
    all_aas = set()
    for i in range(3):
        dna = original_seq[i:]
        dna = dna[:len(dna) - (len(dna) % 3)]
        peptide = translate_dna(dna)
        frames[f'+{i+1}'] = {'dna': dna, 'peptide': peptide}
    
    for i in range(3):
        dna = rev_comp_seq[i:]
        dna = dna[:len(dna) - (len(dna) % 3)]
        peptide = translate_dna(dna)
        frames[f'-{i+1}'] = {'dna': dna, 'peptide': peptide}

    # Gather all amino acids from all frames
    for data in frames.values():
        all_aas.update(list(data['peptide'].replace('*','')))

    print("All Amino Acids found across 6 frames:")
    print(", ".join(sorted(list(all_aas))))
    print("-" * 35)

    # Step 3: Identify the unique frame
    target_frame = None
    target_frame_name = ''
    unique_aas_in_target = set()

    for name, data in frames.items():
        other_peptides = "".join([f['peptide'] for n, f in frames.items() if n != name])
        current_unique = {aa for aa in data['peptide'] if aa not in other_peptides and aa != '*'}
        
        if len(current_unique) == 2:
            aa1, aa2 = list(current_unique)
            pol1 = aa_properties[aa1][1]
            pol2 = aa_properties[aa2][1]
            if (pol1 == 'Polar' and pol2 == 'Non-polar') or (pol1 == 'Non-polar' and pol2 == 'Polar'):
                target_frame = data
                target_frame_name = name
                unique_aas_in_target = current_unique
                break

    print(f"--- Step 2: Identifying the Correct Reading Frame ---")
    if not target_frame:
        print("Error: Could not find a frame that meets the criteria.")
        return
        
    print(f"Found target frame: {target_frame_name}")
    print(f"DNA sequence: {target_frame['dna']}")
    print(f"Peptide sequence: {'-'.join(aa_properties[aa][0] for aa in target_frame['peptide'])}")
    print(f"Unique amino acids: {', '.join(aa_properties[aa][0] for aa in unique_aas_in_target)}")
    
    polar_aa = next(aa for aa in unique_aas_in_target if aa_properties[aa][1] == 'Polar')
    nonpolar_aa = next(aa for aa in unique_aas_in_target if aa_properties[aa][1] == 'Non-polar')
    print(f"This frame contains one unique polar AA ({aa_properties[polar_aa][0]}) and one unique non-polar AA ({aa_properties[nonpolar_aa][0]}), matching the conditions.")
    print("-" * 35)

    # Step 4 & 5: Determine SNPs and construct modified sequence
    original_codons = textwrap.wrap(target_frame['dna'], 3)
    original_peptide = target_frame['peptide']
    
    polar_aa_codon = original_codons[original_peptide.find(polar_aa)]
    nonpolar_aa_codon = original_codons[original_peptide.find(nonpolar_aa)]
    
    # Find new codons after SNP
    modified_polar_codon = next(sc for sc in stop_codons if sum(1 for a, b in zip(polar_aa_codon, sc) if a != b) == 1)
    modified_nonpolar_codon = next(cc for cc in cys_codons if sum(1 for a, b in zip(nonpolar_aa_codon, cc) if a != b) == 1)

    print("--- Step 3: Determining the Specific SNPs ---")
    print(f"Change 1: Polar '{aa_properties[polar_aa][0]}' to STOP codon.")
    print(f"   Original codon: {polar_aa_codon} -> Modified codon: {modified_polar_codon}")
    print(f"Change 2: Non-polar '{aa_properties[nonpolar_aa][0]}' to Cysteine.")
    print(f"   Original codon: {nonpolar_aa_codon} -> Modified codon: {modified_nonpolar_codon}")
    print("-" * 35)
    
    # Step 6: Design the oligo
    modified_codons = [c.replace(polar_aa_codon, modified_polar_codon).replace(nonpolar_aa_codon, modified_nonpolar_codon) for c in original_codons]
    
    stop_index = modified_codons.index(modified_polar_codon)
    oligo_target_seq = "".join(modified_codons[:stop_index])
    final_oligo = get_reverse_complement(oligo_target_seq)
    
    print("--- Step 4: Designing the Oligo ---")
    print("The oligo binds to the modified sequence that is translated into amino acids.")
    print(f"Modified DNA sequence (up to STOP): {' '.join(modified_codons[:stop_index])}")
    print(f"Target sequence for oligo: 5' {oligo_target_seq} 3'")
    print(f"Oligo sequence (reverse complement of target): 5' {final_oligo} 3'")

    print("\nFinal Answer:")
    final_answer = " ".join(textwrap.wrap(final_oligo, 3))
    print(f"<<<5' {final_answer} 3'>>>")

if __name__ == '__main__':
    main()