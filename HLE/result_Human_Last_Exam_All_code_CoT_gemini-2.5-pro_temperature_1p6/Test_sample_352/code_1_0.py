def suggest_mutagenesis_plan():
    """
    This function outlines a plan for a site-directed mutagenesis experiment
    to relieve the negative charge of a specific protein patch.
    """
    print("### Plan for Site-Directed Mutagenesis ###")
    print("Goal: To relieve the autoinhibitory effect of the negatively charged patch at positions 47-50 (S-E-E-D) in protein x.\n")
    print("Strategy: Replace each amino acid in the patch with Alanine (A). Alanine is a small, non-polar, and uncharged amino acid. This 'Alanine scanning' approach is ideal for neutralizing charge and removing functional groups (like the phosphorylation site on Serine) with minimal structural disruption.\n")
    
    print("--- Step-by-Step Rationale ---\n")

    # Original and proposed sequences
    original_positions = [47, 48, 49, 50]
    original_residues = ['S (Serine)', 'E (Glutamate)', 'E (Glutamate)', 'D (Aspartate)']
    proposed_residue = 'A (Alanine)'
    
    # Position 47
    print(f"Position {original_positions[0]}: Original is {original_residues[0]}")
    print(f"  - Problem: The hydroxyl group of Serine is a phosphorylation site, which introduces a strong negative charge.")
    print(f"  - Solution: Replace with {proposed_residue}. Alanine cannot be phosphorylated, thus preventing the addition of a negative charge at this position.")
    print("-" * 30)

    # Position 48
    print(f"Position {original_positions[1]}: Original is {original_residues[1]}")
    print(f"  - Problem: Glutamate is an acidic amino acid with a permanent negative charge.")
    print(f"  - Solution: Replace with {proposed_residue}. This removes the side chain's negative charge.")
    print("-" * 30)
    
    # Position 49
    print(f"Position {original_positions[2]}: Original is {original_residues[2]}")
    print(f"  - Problem: Glutamate is negatively charged.")
    print(f"  - Solution: Replace with {proposed_residue} to continue building a neutral patch.")
    print("-" * 30)

    # Position 50
    print(f"Position {original_positions[3]}: Original is {original_residues[3]}")
    print(f"  - Problem: Aspartate is also an acidic amino acid with a permanent negative charge.")
    print(f"  - Solution: Replace with {proposed_residue} to completely neutralize the patch.")
    print("-" * 30)

    print("\n### Final Recommendation ###")
    print("The best replacement to test the hypothesis that the negative charge of this patch is inhibitory involves mutating all four residues to Alanine.\n")
    
    print("Original sequence equation:")
    print(f"  {original_positions[0]}-{original_residues[0][0]} -> {original_positions[1]}-{original_residues[1][0]} -> {original_positions[2]}-{original_residues[2][0]} -> {original_positions[3]}-{original_residues[3][0]}")

    print("\nProposed new sequence equation:")
    print(f"  {original_positions[0]}-{proposed_residue[0]} -> {original_positions[1]}-{proposed_residue[0]} -> {original_positions[2]}-{proposed_residue[0]} -> {original_positions[3]}-{proposed_residue[0]}")


if __name__ == '__main__':
    suggest_mutagenesis_plan()
<<<AAAA>>>