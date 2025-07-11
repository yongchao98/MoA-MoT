def suggest_mutagenesis():
    """
    Suggests and prints the best amino acid replacements for a site-directed
    mutagenesis experiment based on the provided protein information.
    """
    original_positions = [47, 48, 49, 50]
    original_amino_acids = ["Serine (S)", "Glutamate (E)", "Glutamate (E)", "Aspartate (D)"]

    # Rationale for replacement:
    # The goal is to remove the highly negative charge of the S-E-E-D patch.
    # The standard and most effective approach is to replace charged/functional
    # residues with Alanine (A). Alanine is small, neutral, non-polar, and
    # unlikely to cause steric clashes or introduce new unintended interactions.
    # Replacing Serine-47 with Alanine also crucially removes the phosphorylation site.

    mutant_amino_acid_code = "A"
    mutant_amino_acid_name = "Alanine"

    print("--- Site-Directed Mutagenesis Plan ---")
    print("\nObjective: Relieve the negatively charged environment of the patch at positions 47-50.\n")

    print("Original Patch:")
    original_patch_str = []
    for pos, aa in zip(original_positions, original_amino_acids):
        original_patch_str.append(f"{pos}:{aa.split(' ')[1]}")
    print(f"  Positions 47-50: {' - '.join(original_patch_str)}")
    print("  This patch (pS-E-E-D) is highly negatively charged, causing autoinhibition.\n")


    print("Proposed Mutant Patch:")
    mutant_patch_str = []
    for pos in original_positions:
        mutant_patch_str.append(f"{pos}:{mutant_amino_acid_code}")
    print(f"  Replace with: {' - '.join(mutant_patch_str)}")
    print(f"  The best replacement for each of these four amino acids is Alanine ({mutant_amino_acid_code}).\n")


    print("Final Equation (Original -> Mutant):")
    final_equation_str = ""
    for i in range(len(original_positions)):
        original_aa_code = original_amino_acids[i].split(' ')[1].replace('(', '').replace(')', '')
        final_equation_str += f"{original_positions[i]}: {original_aa_code} -> {mutant_amino_acid_code}"
        if i < len(original_positions) - 1:
            final_equation_str += " | "
            
    print(final_equation_str)


# Execute the function to print the plan
suggest_mutagenesis()
