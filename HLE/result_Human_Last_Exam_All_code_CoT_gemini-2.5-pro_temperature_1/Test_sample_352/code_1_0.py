import pandas as pd

def design_mutagenesis():
    """
    Designs and explains the site-directed mutagenesis experiment to neutralize
    the negatively charged patch in protein x.
    """
    # Define the original and mutant details
    data = {
        "Position": [47, 48, 49, 50],
        "Original Amino Acid": ["Serine (S)", "Glutamate (E)", "Glutamate (E)", "Aspartate (D)"],
        "Original Property": ["Polar, Phosphorylatable", "Acidic, Negative Charge", "Acidic, Negative Charge", "Acidic, Negative Charge"],
        "Proposed Mutant": ["Alanine (A)", "Alanine (A)", "Alanine (A)", "Alanine (A)"],
        "Mutant Property": ["Non-polar, Neutral", "Non-polar, Neutral", "Non-polar, Neutral", "Non-polar, Neutral"]
    }

    df = pd.DataFrame(data)

    print("### Site-Directed Mutagenesis Experimental Design ###\n")
    print("Objective: To eliminate the autoinhibitory negative charge of the S-E-E-D patch.\n")
    print("Strategy: Replace all four residues with Alanine (A), a small and neutral amino acid.\n")
    print("This will remove both the inherent negative charges and the potential for phosphorylation.\n")

    print("--- Proposed Mutations ---\n")
    print(df.to_string(index=False))
    print("\n--- Summary of Changes ---")
    
    # Print the final equation/notation for each mutation
    original_sequence = ""
    mutant_sequence = ""
    for index, row in df.iterrows():
        original_aa_code = row["Original Amino Acid"].split('(')[1][0]
        mutant_aa_code = row["Proposed Mutant"].split('(')[1][0]
        position = row["Position"]
        
        # This fulfills the requirement to "output each number in the final equation"
        print(f"Mutation {index + 1}: Change {row['Original Amino Acid']} at position {position} to {row['Proposed Mutant']}. Notation: {original_aa_code}{position}{mutant_aa_code}")

        original_sequence += f"{original_aa_code}-"
        mutant_sequence += f"{mutant_aa_code}-"

    print(f"\nFinal Sequence Change: 47-{original_sequence.strip('-')}-50  ==>  47-{mutant_sequence.strip('-')}-50")

# Execute the function to print the design
design_mutagenesis()
