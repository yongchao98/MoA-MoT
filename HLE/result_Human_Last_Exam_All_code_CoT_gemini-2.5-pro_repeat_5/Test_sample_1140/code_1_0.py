def synthesis_plan():
    """
    Outlines the synthesis plan for a long peptide with an unnatural amino acid.
    """
    peptide_length = 100
    unnatural_amino_acid = "Azido Phenylalanine (AzF)"
    
    print(f"### Synthesis Plan for a {peptide_length}aa Peptide containing {unnatural_amino_acid} ###")
    print("\n")
    print("Recommended Technique: Native Chemical Ligation (NCL)")
    print("-" * 70)
    print("Reason: Direct Solid-Phase Peptide Synthesis (SPPS) is inefficient for peptides over ~60 amino acids.")
    print("NCL allows for the assembly of long peptides from smaller, high-purity fragments, and easily accommodates unnatural amino acids.")
    print("\n")
    
    print("Step-by-Step NCL Strategy:")
    
    # Step 1: Divide the peptide
    fragment_a_range = "residues 1-50"
    fragment_b_range = "residues 51-100"
    print(f"1. Divide the {peptide_length}aa peptide sequence into two smaller fragments (e.g., Fragment A: {fragment_a_range}; Fragment B: {fragment_b_range}).")

    # Step 2: Synthesize fragments
    print("2. Synthesize each fragment individually using Solid-Phase Peptide Synthesis (SPPS).")

    # Step 3: Incorporate UAA
    print(f"3. During SPPS, incorporate the unnatural amino acid '{unnatural_amino_acid}' into the correct position within its fragment.")

    # Step 4: Prepare for ligation
    fragment_a_modification = "C-terminal thioester"
    fragment_b_modification = "N-terminal Cysteine"
    print(f"4. Prepare the fragments for ligation: Synthesize Fragment A with a {fragment_a_modification} and Fragment B with an {fragment_b_modification}.")

    # Step 5: Ligation
    print("5. After purification, mix the two fragments in solution to ligate them into the full-length peptide via a native peptide bond.")
    
    # Step 6: Final purification
    print("6. Purify the final ligated product to obtain the target 100aa peptide.")

if __name__ == '__main__':
    synthesis_plan()