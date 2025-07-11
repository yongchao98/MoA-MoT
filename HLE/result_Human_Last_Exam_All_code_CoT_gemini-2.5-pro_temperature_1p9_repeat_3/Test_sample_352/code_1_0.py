def design_mutagenesis():
    """
    Designs and prints a site-directed mutagenesis plan to relieve the negative
    charge of a specific patch in protein x.
    """
    # Original amino acid patch: Position, Single-Letter Code, Full Name
    original_patch = [
        (47, 'S', 'Serine'),
        (48, 'E', 'Glutamate'),
        (49, 'E', 'Glutamate'),
        (50, 'D', 'Aspartate')
    ]

    # Proposed replacement amino acid
    replacement_aa = ('A', 'Alanine')

    # Print the experimental plan
    print("--- Site-Directed Mutagenesis Plan for Protein X ---")
    print("Objective: To eliminate the autoinhibitory negative charge from the patch at positions 47-50.\n")
    print("Original Sequence Patch: S-E-E-D (Ser-Glu-Glu-Asp)")
    print("Proposed Mutant Sequence: A-A-A-A (Ala-Ala-Ala-Ala)\n")
    print("--- Individual Mutation Details ---")

    for position, original_code, original_name in original_patch:
        new_code, new_name = replacement_aa
        print(f"Position {position}: Mutate {original_name} ({original_code}) to {new_name} ({new_code})")

    print("\n--- Rationale for Choosing Alanine ---")
    print("1. Charge Neutralization: Alanine is a nonpolar, neutral amino acid. This will remove the negative charges from E48, E49, and D50.")
    print("2. Prevention of Phosphorylation: Replacing Serine at position 47 with Alanine removes the hydroxyl group required for phosphorylation, thus eliminating this potential source of a strong negative charge.")
    print("3. Minimal Structural Disruption: Alanine has a small side chain, minimizing the risk of altering the protein's overall structure, a standard practice known as 'Alanine Scanning'.")

# Execute the function to display the plan
design_mutagenesis()
