def suggest_mutagenesis():
    """
    Suggests amino acid replacements for a site-directed mutagenesis experiment
    to relieve a negatively charged patch in protein x.
    """
    # Define the original amino acid information
    positions = [47, 48, 49, 50]
    original_aa_one_letter = ['S', 'E', 'E', 'D']
    original_aa_three_letter = ['Serine', 'Glutamate', 'Glutamate', 'Aspartate']

    # The ideal replacement to remove charge and phosphorylation potential is Alanine (A)
    replacement_aa = 'A'
    replacement_aa_full_name = 'Alanine'

    print("--- Site-Directed Mutagenesis Plan ---")
    print("Goal: To relieve the autoinhibitory negative charge from the patch at positions 47-50.\n")
    print(f"Original Sequence (47-50): {'-'.join(original_aa_one_letter)}")
    print("This patch is negatively charged and the Serine at position 47 is a phosphorylation site, which adds more negative charge.")
    print(f"Recommendation: Replace all four amino acids with Alanine (A), which is small, neutral, and cannot be phosphorylated.\n")

    print("--- Proposed Mutations ---")
    for i in range(len(positions)):
        original_str = f"{original_aa_three_letter[i]} ({original_aa_one_letter[i]})"
        replacement_str = f"{replacement_aa_full_name} ({replacement_aa})"
        print(f"At position {positions[i]}, replace original {original_str} with new {replacement_str}.")

    final_sequence = replacement_aa * len(positions)
    print("\n--- Final Recommended Sequence (47-50) ---")
    print(final_sequence)

suggest_mutagenesis()
<<<AAAA>>>