def suggest_mutagenesis():
    """
    Suggests and explains the site-directed mutagenesis plan to neutralize
    a negatively charged patch in protein x.
    """
    # Define the original amino acids and their positions
    original_patch = {
        47: "S (Serine)",
        48: "E (Glutamate)",
        49: "E (Glutamate)",
        50: "D (Aspartate)"
    }

    # Define the proposed replacement amino acid
    replacement_aa = "A (Alanine)"
    replacement_sequence = "AAAA"

    # Explain the rationale
    print("Goal: To relieve the autoinhibitory effect of the negatively charged patch at amino acid positions 47-50.\n")
    print("Native Patch Analysis:")
    print("- Positions 48 (E), 49 (E), and 50 (D) are acidic, negatively charged amino acids.")
    print("- Position 47 (S) is a phosphorylation site, which adds a strong negative charge when modified.")
    print("\nProposed Solution:")
    print(f"The best strategy is to replace all four residues with Alanine (A). Alanine is a small, neutral, and non-phosphorylatable amino acid. This 'Alanine scan' will effectively neutralize the entire patch with minimal structural disruption, directly testing the hypothesis that the negative charge is inhibitory.")
    print("\n--- Recommended Mutations ---")

    # Output the mutation for each position
    for position, original_aa in original_patch.items():
        # Extract the single-letter code for the original amino acid
        original_code = original_aa.split(" ")[0]
        # Get the single-letter code for the replacement
        replacement_code = replacement_aa.split(" ")[0]
        print(f"Position {position}: Mutate {original_aa} to {replacement_aa} ({original_code}{position}{replacement_code})")

    print(f"\nThis changes the original sequence 'SEED' at positions 47-50 to '{replacement_sequence}'.")

# Execute the function to print the recommendation
suggest_mutagenesis()