def suggest_mutation():
    """
    Analyzes the protein x amino acid patch and suggests a replacement
    to remove its negative charge for a site-directed mutagenesis experiment.
    """

    # Define the original amino acid patch
    original_patch = {
        47: {'aa': 'Serine', 'code': 'S', 'property': 'Polar, can be phosphorylated (becomes negatively charged)'},
        48: {'aa': 'Glutamate', 'code': 'E', 'property': 'Acidic, negatively charged'},
        49: {'aa': 'Glutamate', 'code': 'E', 'property': 'Acidic, negatively charged'},
        50: {'aa': 'Aspartate', 'code': 'D', 'property': 'Acidic, negatively charged'}
    }

    # Explain the problem
    print("--- Analysis of the Original Patch ---")
    print("The goal is to eliminate the negative charge from the S-E-E-D patch at positions 47-50.")
    print("The negative charge is contributed by:\n"
          "1. Glutamate (E) and Aspartate (D), which are acidic.\n"
          "2. Serine (S), which becomes highly negatively charged upon phosphorylation.\n")

    # Explain the proposed solution
    print("--- Proposed Mutation Strategy ---")
    print("The best strategy is to replace all four residues with Alanine (A).")
    print("Reasoning for choosing Alanine (A):\n"
          "- It is a small and chemically inert amino acid.\n"
          "- Its side chain is neutral, which removes the negative charges from E and D.\n"
          "- It cannot be phosphorylated, which prevents the negative charge addition at position 47.\n"
          "- This 'Alanine scanning' approach minimizes structural disruption while effectively testing the function of the original side chains.\n")
    
    # Define the proposed new amino acid patch
    proposed_patch = {
        47: {'aa': 'Alanine', 'code': 'A'},
        48: {'aa': 'Alanine', 'code': 'A'},
        49: {'aa': 'Alanine', 'code': 'A'},
        50: {'aa': 'Alanine', 'code': 'A'}
    }

    # Construct the final equation string showing the mutation
    original_sequence = "-".join([f"{original_patch[p]['code']}({p})" for p in sorted(original_patch.keys())])
    proposed_sequence = "-".join([f"{proposed_patch[p]['code']}({p})" for p in sorted(proposed_patch.keys())])
    
    print("--- Final Proposed Mutation Equation ---")
    print(f"Original Sequence: {original_sequence}")
    print(f"Proposed Mutant:   {proposed_sequence}")

# Run the function to display the results
suggest_mutation()
