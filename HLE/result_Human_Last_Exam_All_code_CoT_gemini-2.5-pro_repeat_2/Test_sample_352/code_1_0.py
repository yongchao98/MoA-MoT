def propose_mutagenesis():
    """
    Analyzes the S-E-E-D patch and suggests the best amino acid replacements
    to relieve the negative charge for a site-directed mutagenesis experiment.
    """
    # Define the positions and original amino acids of the inhibitory patch
    positions = [47, 48, 49, 50]
    original_sequence = ['S', 'E', 'E', 'D']
    
    # The best replacement amino acid is Alanine ('A') for all positions
    # 'A' is neutral, non-polar, small, and prevents phosphorylation at S47
    mutated_sequence = ['A', 'A', 'A', 'A']
    
    print("Site-Directed Mutagenesis Plan:")
    print("To eliminate the inhibitory negative charge of the S-E-E-D patch,")
    print("the recommended replacement for all four residues is Alanine (A).\n")

    # Build the "equation" strings to clearly show the change
    original_str = " + ".join([f"{aa}{pos}" for pos, aa in zip(positions, original_sequence)])
    mutated_str = " + ".join([f"{aa}{pos}" for pos, aa in zip(positions, mutated_sequence)])

    print("Proposed Mutational Change:")
    print(f"{original_str}  ->  {mutated_str}")

# Execute the function to print the proposal
propose_mutagenesis()
