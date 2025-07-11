def suggest_mutagenesis():
    """
    Prints the recommended amino acid substitutions for the site-directed
    mutagenesis experiment on protein x.
    """
    # Original amino acids and their positions
    original_patch = {
        47: 'S (Serine)',
        48: 'E (Glutamate)',
        49: 'E (Glutamate)',
        50: 'D (Aspartate)'
    }

    # Proposed replacement amino acid
    replacement_aa = 'A (Alanine)'
    
    print("--- Site-Directed Mutagenesis Plan ---")
    print("Objective: To relieve the autoinhibitory effect of the negatively charged patch (aa 47-50).\n")
    print("Strategy: Replace the negatively charged and phosphorylatable residues with a small, neutral amino acid to eliminate the negative charge while minimizing structural disruption.\n")
    print("Recommendation: Mutate all four residues to Alanine (A).\n")
    print("--- Proposed Mutation Details ---")

    # Print the mapping from original to proposed amino acids
    for position, aa_name in original_patch.items():
        print(f"Position {position}: {aa_name} -> {replacement_aa}")

    print("\n--- Final Recommended Sequence ---")
    final_sequence = "47-A-A-A-A-50"
    print(f"The proposed mutant sequence for positions 47-50 is: {final_sequence}")

if __name__ == '__main__':
    suggest_mutagenesis()