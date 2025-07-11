def suggest_mutagenesis():
    """
    Provides a site-directed mutagenesis plan to relieve a negatively charged
    inhibitory patch in protein x.
    """
    # Define the positions and the original amino acid sequence segment
    positions = [47, 48, 49, 50]
    original_amino_acids = ['S', 'E', 'E', 'D']
    original_sequence = "".join(original_amino_acids)

    # Explain the rationale for the mutation
    print("### Site-Directed Mutagenesis Plan for Protein X ###\n")
    print(f"Objective: To neutralize the autoinhibitory patch at positions {min(positions)}-{max(positions)}.")
    print(f"Original Sequence Segment: {original_sequence}\n")

    print("Rationale for replacements:")
    print("- Position 47 (Serine, S): A phosphorylation site. When phosphorylated, it adds a strong negative charge.")
    print("- Positions 48, 49 (Glutamate, E) and 50 (Aspartate, D): These are acidic amino acids, each carrying a negative charge.")
    print("\nProposed Solution: Replace all four residues with Alanine (A).")
    print("- Alanine is small, neutral, and cannot be phosphorylated.")
    print("- This mutation will effectively remove the negative charges and prevent phosphorylation.\n")

    # Define the proposed mutant amino acids
    mutant_amino_acids = ['A', 'A', 'A', 'A']
    mutant_sequence = "".join(mutant_amino_acids)

    # Display the specific changes for the experiment
    print("--- Experimental Changes ---")
    for i in range(len(positions)):
        print(f"Position {positions[i]}: Change original '{original_amino_acids[i]}' to proposed '{mutant_amino_acids[i]}'")
    print("--------------------------\n")
    print(f"Final Result: The original sequence '{original_sequence}' should be mutated to '{mutant_sequence}'.")


if __name__ == '__main__':
    suggest_mutagenesis()
