def find_oligo_sequence():
    """
    Solves the DNA puzzle by identifying the correct reading frame, applying the specified SNPs,
    and designing the complementary oligo for the resulting translated sequence.
    """

    # --- Step 1: Define the sequence and codon/amino acid information ---
    original_sequence = "CTT CCC CGC ACA AGT GGT"
    
    # In our analysis, we identified Frame +2 as the correct frame.
    # The original sequence is 5' CTT CCC CGC ACA AGT GGT 3'
    # Frame +2 starts from the second nucleotide.
    frame_plus_2_dna = "TTC CCG CAC AAG TGG"
    
    codon_table = {
        'TTC': 'F', 'CCG': 'P', 'CAC': 'H', 'AAG': 'K', 'TGG': 'W',
        'TGC': 'C', 'TAG': 'Stop'
    }
    
    amino_acid_properties = {
        'F': 'Non-polar', 'P': 'Non-polar', 'H': 'Polar', 
        'K': 'Polar', 'W': 'Non-polar', 'C': 'Cysteine'
    }

    print("--- Step-by-step analysis ---")
    print(f"Original Sequence: 5' {original_sequence} 3'")
    print(f"Target Reading Frame (Frame +2): 5' {frame_plus_2_dna} 3'")

    # --- Step 2: Explain the SNP selection based on the problem's rules ---
    original_aa_sequence = " ".join([codon_table[c] for c in frame_plus_2_dna.split()])
    print(f"Original Amino Acid Translation: {original_aa_sequence}\n")
    
    print("Identifying the two required Single Nucleotide Polymorphisms (SNPs):")
    # SNP 1: Polar amino acid to a Stop codon
    polar_codon_original = "AAG"
    polar_aa_original = codon_table[polar_codon_original]
    polar_codon_modified = "TAG"
    stop_signal = codon_table[polar_codon_modified]
    print(f"1. A polar amino acid is changed to a STOP codon.")
    print(f"   - The codon '{polar_codon_original}' (for {polar_aa_original}, {amino_acid_properties[polar_aa_original]}) is changed to '{polar_codon_modified}' ({stop_signal}).")

    # SNP 2: Non-polar amino acid to a Cysteine amino acid
    nonpolar_codon_original = "TTC"
    nonpolar_aa_original = codon_table[nonpolar_codon_original]
    nonpolar_codon_modified = "TGC"
    cysteine_aa = codon_table[nonpolar_codon_modified]
    print(f"2. A non-polar amino acid is changed to a Cysteine.")
    print(f"   - The codon '{nonpolar_codon_original}' (for {nonpolar_aa_original}, {amino_acid_properties[nonpolar_aa_original]}) is changed to '{nonpolar_codon_modified}' ({cysteine_aa}).\n")

    # --- Step 3: Determine the sequence to be targeted by the oligo ---
    # The oligo binds to the part of the modified sequence that is translated into amino acids.
    # This means the sequence before the new stop codon.
    # Modified DNA: TGC CCG CAC TAG TGG
    target_sequence_str = f"{nonpolar_codon_modified} CCG CAC"
    modified_codons = target_sequence_str.split()
    translated_aas = " ".join([codon_table[c] for c in modified_codons])
    
    print("Constructing the target sequence for the oligo:")
    print(f"The modified DNA sequence segment that gets translated is 5' {target_sequence_str} 3'")
    print(f"This translates to the peptide: {translated_aas}\n")

    # --- Step 4: Calculate the reverse complement to get the oligo sequence ---
    target_dna = target_sequence_str.replace(" ", "")
    
    complement_map = str.maketrans('ATCG', 'TAGC')
    complement_dna = target_dna.translate(complement_map)
    oligo_sequence = complement_dna[::-1]
    
    print("Designing the oligo:")
    print("The oligo sequence must be the reverse complement of the target sequence.")
    print(f"Target:         5' {target_dna} 3'")
    print(f"Complement:     3' {complement_dna} 5'")
    print(f"Final Oligo (Reverse Complement): 5' {oligo_sequence} 3'")

    return oligo_sequence

# --- Final Answer ---
final_answer = find_oligo_sequence()
print(f"\n<<<5' {final_answer} 3'>>>")
