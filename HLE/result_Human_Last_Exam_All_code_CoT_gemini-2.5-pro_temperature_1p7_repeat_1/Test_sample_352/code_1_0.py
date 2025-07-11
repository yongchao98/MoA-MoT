import pandas as pd

def propose_mutagenesis():
    """
    Proposes amino acid substitutions to neutralize a negatively charged patch in a protein.
    """
    # Define the original amino acids and their properties
    original_patch = {
        'Position': [47, 48, 49, 50],
        'Original Amino Acid': ['Serine (S)', 'Glutamate (E)', 'Glutamate (E)', 'Aspartate (D)'],
        'Property': ['Phosphorylation site, becomes negatively charged', 'Negatively charged', 'Negatively charged', 'Negatively charged']
    }

    # The chosen replacement is Alanine (A) for its neutral and small side chain.
    replacement_amino_acid = 'Alanine (A)'

    print("Goal: To relieve the autoinhibitory negative charge from the patch at positions 47-50.")
    print("Strategy: Replace the native residues with a small, neutral amino acid to eliminate charge with minimal structural disruption.")
    print(f"The best choice for this is {replacement_amino_acid}.\n")

    print("--- Proposed Mutations ---")
    
    # Create a DataFrame for clear presentation
    df = pd.DataFrame(original_patch)
    df['Proposed Replacement'] = replacement_amino_acid
    
    # Print each mutation individually as requested
    for index, row in df.iterrows():
        print(f"At Position {row['Position']}: Replace {row['Original Amino Acid']} with {row['Proposed Replacement']}")

    print("\nThe original sequence S-E-E-D should be mutated to A-A-A-A.")

# Run the function to display the plan and proposed mutations
propose_mutagenesis()
