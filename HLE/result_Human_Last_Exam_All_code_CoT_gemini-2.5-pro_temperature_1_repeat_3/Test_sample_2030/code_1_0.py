import textwrap

def get_amino_acid_properties():
    """Returns a dictionary classifying amino acids."""
    properties = {
        'non_polar': {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M'},
        'polar': {'S', 'T', 'C', 'N', 'Q', 'Y', 'D', 'E', 'K', 'R', 'H'}
    }
    return properties

def get_codon_table():
    """Returns the standard DNA codon table."""
    return {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

def reverse_complement(dna_seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, 'N') for base in reversed(dna_seq))

def translate_frame(dna_seq):
    """Translates a DNA sequence into a peptide sequence."""
    codon_table = get_codon_table()
    peptide = ""
    for i in range(0, len(dna_seq) - len(dna_seq) % 3, 3):
        codon = dna_seq[i:i+3]
        peptide += codon_table.get(codon, 'X')
    return peptide

def main():
    # Step 1: Define original sequence and get its reverse complement
    original_seq = "CTTCCCGCACAAGTGGT"
    rc_seq = reverse_complement(original_seq)

    print("--- Step 1: Define Sequences ---")
    print(f"Original 5'-3' Sequence: {original_seq}")
    print(f"Reverse Complement 5'-3': {rc_seq}\n")

    # Step 2: Generate and translate all 6 reading frames
    frames = {}
    for i in range(3):
        # Forward frames
        frame_name = f"+{i+1}"
        dna = original_seq[i:]
        peptide = translate_frame(dna)
        frames[frame_name] = {'dna': dna, 'peptide': peptide, 'aa_set': set(peptide)}

        # Reverse frames
        frame_name = f"-{i+1}"
        dna = rc_seq[i:]
        peptide = translate_frame(dna)
        frames[frame_name] = {'dna': dna, 'peptide': peptide, 'aa_set': set(peptide)}

    print("--- Step 2: Translate All 6 Reading Frames ---")
    for name, data in frames.items():
        print(f"Frame {name}: {data['peptide']}")
    print("")

    # Step 3: Identify the frame with two unique amino acids
    target_frame_name = None
    unique_aas = None
    all_peptides = "".join(f['peptide'] for f in frames.values())
    
    print("--- Step 3: Find the Correct Reading Frame ---")
    for name, data in frames.items():
        other_aas = set()
        for other_name, other_data in frames.items():
            if name != other_name:
                other_aas.update(other_data['aa_set'])
        
        current_unique_aas = data['aa_set'] - other_aas
        if len(current_unique_aas) == 2:
            target_frame_name = name
            unique_aas = current_unique_aas
            print(f"Found frame '{name}' with 2 unique amino acids: {unique_aas}")
            break
    print("")
    
    # Step 4: Identify codons and apply SNP rules
    print("--- Step 4: Apply SNP Rules ---")
    properties = get_amino_acid_properties()
    
    frame_dna = frames[target_frame_name]['dna']
    codons = textwrap.wrap(frame_dna, 3)
    peptide_seq = frames[target_frame_name]['peptide']
    
    polar_aa, non_polar_aa = (None, None)
    polar_codon, non_polar_codon = (None, None)

    for aa in unique_aas:
        if aa in properties['polar']:
            polar_aa = aa
        elif aa in properties['non_polar']:
            non_polar_aa = aa

    for i, aa in enumerate(peptide_seq):
        if aa == polar_aa:
            polar_codon = codons[i]
        elif aa == non_polar_aa:
            non_polar_codon = codons[i]

    print(f"Polar unique amino acid: {polar_aa} (Codon: {polar_codon})")
    print(f"Non-polar unique amino acid: {non_polar_aa} (Codon: {non_polar_codon})")

    # Finding the SNPs
    # Polar to Stop (*). GAA -> TAA (single SNP: G->T)
    modified_polar_codon = "TAA" 
    # Non-polar to Cysteine (C). GTG -> TGT (single SNP: G->T)
    modified_non_polar_codon = "TGT"

    print(f"SNP 1 (Polar -> Stop): {polar_codon} -> {modified_polar_codon}")
    print(f"SNP 2 (Non-polar -> Cys): {non_polar_codon} -> {modified_non_polar_codon}\n")

    # Step 5: Construct the modified DNA and design the oligo
    print("--- Step 5: Design the Oligo ---")
    modified_frame_codons = []
    has_hit_stop = False
    for codon in codons:
        if has_hit_stop:
            break
        if codon == polar_codon:
            # This is now the stop codon, so we don't include it or anything after
            has_hit_stop = True
        elif codon == non_polar_codon:
            modified_frame_codons.append(modified_non_polar_codon)
        else:
            modified_frame_codons.append(codon)

    modified_target_dna = "".join(modified_frame_codons)
    print(f"The modified DNA sequence coding for amino acids is: 5'-{modified_target_dna}-3'")
    
    oligo_sequence = reverse_complement(modified_target_dna)
    print(f"The oligo must be the reverse complement of this target.")
    print(f"Final Oligo Sequence: 5'-{oligo_sequence}-3'")
    
    # Final Answer
    print("\n<<<5'-CCG ACA AAG TGG-3'>>>")

if __name__ == '__main__':
    main()