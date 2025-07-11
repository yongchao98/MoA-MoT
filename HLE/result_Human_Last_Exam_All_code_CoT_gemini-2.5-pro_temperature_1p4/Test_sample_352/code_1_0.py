import pandas as pd

def design_mutagenesis():
    """
    Designs and explains the site-directed mutagenesis experiment
    to relieve the inhibitory, negatively charged patch in protein x.
    """

    # 1. Define the original and replacement amino acids
    original_patch = {
        'Position': [47, 48, 49, 50],
        'AA_3_Letter': ['Ser', 'Glu', 'Glu', 'Asp'],
        'AA_1_Letter': ['S', 'E', 'E', 'D'],
        'Property': [
            'Neutral, but can be phosphorylated to become negative',
            'Negatively charged (acidic)',
            'Negatively charged (acidic)',
            'Negatively charged (acidic)'
        ]
    }

    replacement_patch = {
        'Position': [47, 48, 49, 50],
        'AA_3_Letter': ['Ala', 'Ala', 'Ala', 'Ala'],
        'AA_1_Letter': ['A', 'A', 'A', 'A'],
        'Property': [
            'Neutral, cannot be phosphorylated',
            'Neutral, small side chain',
            'Neutral, small side chain',
            'Neutral, small side chain'
        ]
    }

    # Create dataframes for clear display
    original_df = pd.DataFrame(original_patch)
    replacement_df = pd.DataFrame(replacement_patch)

    # 2. Print the rationale
    print("--- Site-Directed Mutagenesis Experiment Design ---")
    print("\nObjective: Relieve the autoinhibitory effect of the negatively charged patch (aa 47-50).")
    print("\nHypothesis: The negative charges from S(p), E, E, and D are responsible for inhibition.")
    print("\nStrategy: Replace all four amino acids with Alanine (A), a small and neutral amino acid.")
    print("\nRationale for Alanine replacement:")
    print("  - S47A: Prevents phosphorylation, removing a key potential negative charge.")
    print("  - E48A, E49A, D50A: Replaces inherently negatively charged residues with neutral ones.")
    print("  - Overall: This creates a 'neutral linker' to test the charge hypothesis with minimal other structural or chemical changes.")

    # 3. Print the final equation/summary of changes
    print("\n--- Summary of Proposed Mutations ---")
    for i in range(len(original_patch['Position'])):
        pos = original_patch['Position'][i]
        orig_aa = original_patch['AA_1_Letter'][i]
        rep_aa = replacement_patch['AA_1_Letter'][i]
        print(f"Position {pos}: {orig_aa} -> {rep_aa}")

    # Combine into a final "equation" format
    original_sequence = "".join(original_patch['AA_1_Letter'])
    mutated_sequence = "".join(replacement_patch['AA_1_Letter'])
    print(f"\nFinal Proposed Change: {original_sequence} -> {mutated_sequence}")

# Execute the function to print the design
design_mutagenesis()
