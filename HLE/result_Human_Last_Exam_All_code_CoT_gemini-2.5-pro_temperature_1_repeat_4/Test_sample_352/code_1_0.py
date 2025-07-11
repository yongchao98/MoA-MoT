def design_mutagenesis_experiment():
    """
    This script outlines the rationale and design for a site-directed mutagenesis
    experiment to eliminate the inhibitory effect of a negatively charged
    amino acid patch in protein x.
    """
    
    # Define the original state of the amino acid patch
    positions = [47, 48, 49, 50]
    original_amino_acids = {
        47: "Serine (S)",
        48: "Glutamate (E)",
        49: "Glutamate (E)",
        50: "Aspartate (D)"
    }
    original_sequence = "S-E-E-D"

    # Define the proposed changes and the rationale
    # The best replacement is Alanine (A) for all positions.
    proposed_amino_acid = "Alanine (A)"
    proposed_sequence = "A-A-A-A"

    print("--- Mutagenesis Plan to Relieve Protein X Autoinhibition ---")
    print("\nObjective: To eliminate the negative charge from the S-E-E-D patch at positions 47-50.")
    print(f"This will test the hypothesis that this patch is autoinhibitory.\n")
    
    print(f"Original sequence at positions {positions[0]}-{positions[-1]}: {original_sequence}\n")

    print("Rationale for mutations:")
    
    # Position 47
    print(f"1. Position {positions[0]} ({original_amino_acids[positions[0]]}):")
    print("   - This Serine is a phosphorylation site. Phosphorylation adds a strong negative charge.")
    print(f"   - Recommendation: Mutate Serine to Alanine (S{positions[0]}A).")
    print(f"   - Justification: Alanine lacks the hydroxyl group necessary for phosphorylation, thus preventing the addition of a negative charge at this site.\n")

    # Positions 48, 49, 50
    print(f"2. Positions {positions[1]}, {positions[2]}, {positions[3]} ({original_amino_acids[positions[1]]}, {original_amino_acids[positions[2]]}, {original_amino_acids[positions[3]]}):")
    print("   - These are acidic amino acids, which are negatively charged at physiological pH.")
    print(f"   - Recommendation: Mutate all three residues to Alanine (E{positions[1]}A, E{positions[2]}A, D{positions[3]}A).")
    print("   - Justification: Alanine is a small, neutral amino acid. Replacing the charged residues with Alanine will completely neutralize the patch. This 'Alanine scanning' approach is a robust method to ablate the function of specific residues.\n")

    print("--- Summary of Proposed Changes ---")
    print(f"Original Patch ({positions[0]}-{positions[-1]}): {original_amino_acids[positions[0]]} - {original_amino_acids[positions[1]]} - {original_amino_acids[positions[2]]} - {original_amino_acids[positions[3]]}")
    print(f"The final proposed mutant sequence is: {proposed_amino_acid} - {proposed_amino_acid} - {proposed_amino_acid} - {proposed_amino_acid}")
    print("\nFinal Equation (Original -> Mutant):")
    print(f"Position {positions[0]}: {original_sequence[0]} -> {proposed_sequence[0]}")
    print(f"Position {positions[1]}: {original_sequence[2]} -> {proposed_sequence[2]}")
    print(f"Position {positions[2]}: {original_sequence[4]} -> {proposed_sequence[4]}")
    print(f"Position {positions[3]}: {original_sequence[6]} -> {proposed_sequence[6]}")


# Execute the function to print the plan
design_mutagenesis_experiment()