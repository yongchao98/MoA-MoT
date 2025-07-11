def suggest_mutations():
    """
    Suggests mutations to neutralize a negatively charged protein patch.
    """
    # Define the positions and original amino acids
    positions = [47, 48, 49, 50]
    original_aa_info = {
        47: "Serine (S)",
        48: "Glutamate (E)",
        49: "Glutamate (E)",
        50: "Aspartate (D)"
    }
    
    # Define the proposed mutant amino acid
    mutant_aa_info = {
        47: "Alanine (A)",
        48: "Alanine (A)",
        49: "Alanine (A)",
        50: "Alanine (A)"
    }
    
    # Explain the rationale
    print("--- Site-Directed Mutagenesis Plan ---")
    print("Goal: To relieve the autoinhibitory negative charge at the 47-50 patch of protein x.\n")
    print("Strategy: Replace the phosphorylatable and acidic residues with Alanine (A). Alanine is small, neutral, and removes the specific side-chain functionalities responsible for the negative charge.\n")
    
    # Print the proposed changes step-by-step
    print("Specific Mutations:")
    for pos in positions:
        print(f"Position {pos}: Replace original {original_aa_info[pos]} with mutant {mutant_aa_info[pos]}.")

    # Summarize the final result in an "equation" format
    original_sequence = "".join([info.split('(')[1][0] for pos, info in sorted(original_aa_info.items())])
    mutant_sequence = "".join([info.split('(')[1][0] for pos, info in sorted(mutant_aa_info.items())])
    
    print("\n--- Summary ---")
    print(f"To test the hypothesis, mutate the original sequence {original_sequence} at positions {positions[0]}-{positions[-1]} to the neutral sequence {mutant_sequence}.")
    
suggest_mutations()
<<<AAAA>>>