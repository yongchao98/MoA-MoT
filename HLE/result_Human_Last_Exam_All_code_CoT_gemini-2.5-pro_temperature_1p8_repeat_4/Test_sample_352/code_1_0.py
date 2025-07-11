import pandas as pd

def suggest_mutagenesis():
    """
    Suggests and prints the site-directed mutagenesis plan to neutralize
    the negatively charged patch in protein x.
    """

    # Define the original and proposed mutations
    mutations = {
        'Position': [47, 48, 49, 50],
        'Original AA': ['Serine', 'Glutamate', 'Glutamate', 'Aspartate'],
        'Original Code': ['S', 'E', 'E', 'D'],
        'Reason for Negative Charge': [
            'Phosphorylation Site (-PO4)',
            'Acidic Side Chain (-COO)',
            'Acidic Side Chain (-COO)',
            'Acidic Side Chain (-COO)'
        ],
        'Proposed AA': ['Alanine', 'Alanine', 'Alanine', 'Alanine'],
        'Proposed Code': ['A', 'A', 'A', 'A'],
        'Justification': [
            'Cannot be phosphorylated; Neutral',
            'Neutral and small',
            'Neutral and small',
            'Neutral and small'
        ]
    }

    # Create a DataFrame for clear output
    df = pd.DataFrame(mutations)

    print("Experimental Design: Site-Directed Mutagenesis of the S-E-E-D patch")
    print("-" * 70)
    print("Goal: To relieve the autoinhibitory effect by removing the concentrated negative charge.")
    print("Strategy: Replace the four-amino acid patch with Alanine (A) to create a neutral region.\n")

    print(df.to_string(index=False))
    
    print("\nSummary of Proposed Mutation:")
    original_patch = '-'.join(mutations['Original Code'])
    mutant_patch = '-'.join(mutations['Proposed Code'])
    print(f"Original sequence (47-50): {original_patch}")
    print(f"Mutant sequence (47-50):   {mutant_patch}")
    
    final_answer = "The best replacement is to mutate the entire patch S-E-E-D at positions 47-50 to A-A-A-A (Alanine-Alanine-Alanine-Alanine)."
    print(f"\nFinal Recommendation: {final_answer}")


if __name__ == '__main__':
    suggest_mutagenesis()
    # The final answer in the required format
    print("\n<<<The best replacement is to mutate the S-E-E-D patch (positions 47-50) to A-A-A-A (Alanine-Alanine-Alanine-Alanine).>>>")