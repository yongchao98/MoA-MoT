import pandas as pd

def suggest_mutagenesis():
    """
    Analyzes the S-E-E-D patch and suggests the best replacement amino acids
    for a site-directed mutagenesis experiment.
    """
    # Define the original amino acids, their positions, and properties
    original_patch = {
        'Position': [47, 48, 49, 50],
        'Original_AA': ['S', 'E', 'E', 'D'],
        'Name': ['Serine', 'Glutamate', 'Glutamate', 'Aspartate'],
        'Property': ['Polar, Phosphorylatable', 'Negatively Charged', 'Negatively Charged', 'Negatively Charged']
    }

    # The best replacement is Alanine (A) for all positions
    replacement_aa = 'A'
    replacement_name = 'Alanine'
    replacement_property = 'Non-polar, Neutral, Small'

    print("### Site-Directed Mutagenesis Experimental Design ###\n")
    print(f"Objective: To relieve the autoinhibitory effect of the negatively charged patch at positions 47-50 (S-E-E-D).\n")
    print(f"Strategy: Replace all four residues with Alanine (A). Alanine is small, uncharged, and cannot be phosphorylated, making it the ideal choice to neutralize the patch while minimizing structural disruption.\n")

    print("--- Recommended Mutations ---")
    for i in range(len(original_patch['Position'])):
        pos = original_patch['Position'][i]
        orig_aa = original_patch['Original_AA'][i]
        
        # This line constructs and prints the final mutation equation for each position.
        # It shows the position number, original amino acid, and the proposed replacement.
        print(f"Mutation for position {pos}: Change {original_patch['Name'][i]} ({orig_aa}) to {replacement_name} ({replacement_aa}).  Notation: {orig_aa}{pos}{replacement_aa}")

    print("\n--- Summary ---")
    original_sequence = "-".join(original_patch['Original_AA'])
    mutant_sequence = "-".join([replacement_aa] * 4)
    print(f"This experiment will change the amino acid sequence from '{original_sequence}' to '{mutant_sequence}' at positions 47-50.")
    print("This new patch will be chemically neutral and non-phosphorylatable, directly testing the hypothesis of charge-based inhibition.")

# Execute the function to print the recommendation
suggest_mutagenesis()