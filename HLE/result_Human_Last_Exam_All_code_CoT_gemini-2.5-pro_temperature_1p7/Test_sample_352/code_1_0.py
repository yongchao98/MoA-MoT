def suggest_mutations():
    """
    Suggests mutations for a target protein patch to relieve a negative charge.
    """

    # Define the target patch: position, original amino acid, and its properties
    original_patch = {
        47: {'aa': 'Serine', 'code': 'S', 'property': 'Polar, can be phosphorylated to become highly negative'},
        48: {'aa': 'Glutamate', 'code': 'E', 'property': 'Acidic, negatively charged'},
        49: {'aa': 'Glutamate', 'code': 'E', 'property': 'Acidic, negatively charged'},
        50: {'aa': 'Aspartate', 'code': 'D', 'property': 'Acidic, negatively charged'}
    }

    # Define the proposed replacement
    replacement_aa = {'aa': 'Alanine', 'code': 'A', 'property': 'Neutral, non-polar, cannot be phosphorylated'}

    # --- Print the report ---
    print("Plan for Site-Directed Mutagenesis of Protein X (Positions 47-50)")
    print("=" * 70)
    print("Goal: To relieve the autoinhibitory negative charge of the S-E-E-D patch.\n")
    print("Strategy: Replace all four residues with Alanine (A).")
    print(f"Reasoning: {replacement_aa['aa']} ({replacement_aa['code']}) is ideal because it is small, neutral, and")
    print("cannot be phosphorylated. This systematically removes all negative charges from the patch.\n")

    print("--- Proposed Mutations ---")
    mutant_sequence = []
    original_sequence = []
    for position, info in original_patch.items():
        original_sequence.append(info['code'])
        mutant_sequence.append(replacement_aa['code'])
        print(
            f"Position {position}: Replace original {info['aa']} ({info['code']}) "
            f"with {replacement_aa['aa']} ({replacement_aa['code']}). "
            f"Mutation: {info['code']}{position}{replacement_aa['code']}"
        )
    
    print("\n--- Final Sequence Equation ---")
    original_str = "".join(original_sequence)
    mutant_str = "".join(mutant_sequence)
    print(f"The original sequence patch '{original_str}' at positions 47-50 should be mutated to '{mutant_str}'.")
    print(f"Equation: {original_sequence[0]}47 {original_sequence[1]}48 {original_sequence[2]}49 {original_sequence[3]}50 -> {mutant_sequence[0]}47 {mutant_sequence[1]}48 {mutant_sequence[2]}49 {mutant_sequence[3]}50")

suggest_mutations()
<<<AAAA>>>