def suggest_mutagenesis():
    """
    Suggests and explains the optimal amino acid substitutions to neutralize
    the negatively charged S-E-E-D patch in protein x.
    """
    # Original patch information
    positions = [47, 48, 49, 50]
    original_aa = ['S', 'E', 'E', 'D']
    original_names = ['Serine', 'Glutamate', 'Glutamate', 'Aspartate']
    original_properties = [
        "Neutral, but can be phosphorylated to become strongly negative",
        "Negatively charged (acidic)",
        "Negatively charged (acidic)",
        "Negatively charged (acidic)"
    ]

    # Proposed mutant patch information
    mutant_aa = ['A', 'Q', 'Q', 'N']
    mutant_names = ['Alanine', 'Glutamine', 'Glutamine', 'Asparagine']
    mutant_properties = [
        "Neutral, non-phosphorylatable",
        "Neutral, structurally similar to Glutamate",
        "Neutral, structurally similar to Glutamate",
        "Neutral, structurally similar to Aspartate"
    ]

    print("--- Plan for Site-Directed Mutagenesis ---")
    print("Goal: To neutralize the negatively charged patch at positions 47-50.\n")

    print("Original Patch:")
    original_sequence_str = "-".join(original_aa)
    print(f"  Sequence: {original_sequence_str}")
    for i in range(len(positions)):
        print(f"  Position {positions[i]}: {original_names[i]} ({original_aa[i]}) - {original_properties[i]}")

    print("\n-------------------------------------------------\n")

    print("Proposed Mutant Patch:")
    mutant_sequence_str = "-".join(mutant_aa)
    print(f"  Sequence: {mutant_sequence_str}")
    for i in range(len(positions)):
        print(f"  Position {positions[i]}: Replace {original_names[i]} ({original_aa[i]}) with {mutant_names[i]} ({mutant_aa[i]})")
        print(f"  Rationale: Results in a {mutant_properties[i]}.")

    print("\n--- Final Recommendation ---")
    print("The original sequence S-E-E-D at positions 47-50 should be mutated to A-Q-Q-N.")
    print("This change will create a neutral patch to effectively test the hypothesis that the negative charge is responsible for autoinhibition.")

suggest_mutagenesis()
<<<A-Q-Q-N>>>