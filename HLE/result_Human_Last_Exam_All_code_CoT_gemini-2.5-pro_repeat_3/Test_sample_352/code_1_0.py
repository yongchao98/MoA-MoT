def suggest_mutagenesis():
    """
    Suggests and explains the best amino acid replacements for a site-directed
    mutagenesis experiment to neutralize a negatively charged protein patch.
    """

    # --- Define the Original Patch ---
    original_positions = [47, 48, 49, 50]
    original_sequence = ['S', 'E', 'E', 'D']
    original_properties = {
        'S': 'Polar, Phosphorylation Site (becomes -2 charge)',
        'E': 'Acidic, Negatively Charged (-1)',
        'D': 'Acidic, Negatively Charged (-1)'
    }

    # --- Define the Proposed Mutant Patch ---
    mutant_sequence = ['A', 'A', 'A', 'A']
    mutant_properties = {
        'A': 'Non-polar, Neutral Charge'
    }

    # --- Print the Rationale ---
    print("--- Site-Directed Mutagenesis Plan ---")
    print("\nObjective: Relieve the inhibitory effect of the negatively charged patch at positions 47-50.\n")

    # 1. Analyze the original sequence
    print("1. Original Amino Acid Patch:")
    original_str = ""
    for i in range(len(original_positions)):
        pos = original_positions[i]
        aa = original_sequence[i]
        prop = original_properties.get(aa, "Unknown")
        original_str += f"{aa}{pos}"
        if i < len(original_positions) - 1:
            original_str += "-"
        print(f"  - Position {pos}: {aa} ({prop})")
    
    print(f"\nThis original patch '{original_str}' is highly acidic and negatively charged,")
    print("especially when Serine-47 is phosphorylated, which contributes to autoinhibition.\n")

    # 2. Propose the new sequence
    print("2. Proposed Mutant Patch:")
    mutant_str = ""
    for i in range(len(original_positions)):
        pos = original_positions[i]
        aa = mutant_sequence[i]
        prop = mutant_properties.get(aa, "Unknown")
        mutant_str += f"{aa}{pos}"
        if i < len(original_positions) - 1:
            mutant_str += "-"
        print(f"  - Position {pos}: Replace with {aa} ({prop})")

    print("\nRationale for choosing Alanine (A):")
    print("  - It is small, non-polar, and has a neutral charge.")
    print("  - It directly neutralizes the negative charges from Glutamate (E) and Aspartate (D).")
    print("  - Replacing Serine (S) with Alanine (A) prevents phosphorylation at position 47.")
    print("  - This 'Alanine Scan' minimizes structural disruption while testing the function of the charged side chains.\n")
    
    # 3. Final recommendation
    print("--- Final Recommended Mutation ---")
    final_mutation_str = f"Replace original sequence {original_sequence[0]}{original_positions[0]}-{original_sequence[1]}{original_positions[1]}-{original_sequence[2]}{original_positions[2]}-{original_sequence[3]}{original_positions[3]} with {mutant_sequence[0]}{original_positions[0]}-{mutant_sequence[1]}{original_positions[1]}-{mutant_sequence[2]}{original_positions[2]}-{mutant_sequence[3]}{original_positions[3]}"
    print("The recommended mutations are:")
    for i in range(len(original_positions)):
        print(f"  {original_sequence[i]}{original_positions[i]}A (Replace {original_sequence[i]} at position {original_positions[i]} with Alanine)")

# Execute the function
suggest_mutagenesis()
