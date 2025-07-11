import pandas as pd
from io import StringIO

def design_mutagenesis_experiment():
    """
    Designs and explains a site-directed mutagenesis experiment to neutralize
    a negatively charged patch in a protein.
    """

    print("Step 1: Define the target region and the problem.")
    print("The target is the amino acid patch at positions 47-50.")
    print("Original sequence: Serine-Glutamate-Glutamate-Aspartate (S-E-E-D).")
    print("Problem: This region is highly negatively charged, especially when S47 is phosphorylated, leading to suspected autoinhibition.\n")

    print("Step 2: Propose a solution based on the hypothesis.")
    print("The best way to test if the negative charge is inhibitory is to replace")
    print("the charged/phosphorylatable amino acids with a small, neutral amino acid.")
    print("The ideal candidate for this is Alanine (A).\n")

    print("Step 3: Detail the proposed mutations.")
    # Data for the table
    mutation_data = {
        'Position': [47, 48, 49, 50],
        'Original AA (Code)': ['S', 'E', 'E', 'D'],
        'Original AA (Name)': ['Serine', 'Glutamate', 'Glutamate', 'Aspartate'],
        'Reason for Mutation': ['Phosphorylation site, becomes negative', 'Negatively charged', 'Negatively charged', 'Negatively charged'],
        'Proposed Mutant (AA)': ['A (Alanine)', 'A (Alanine)', 'A (Alanine)', 'A (Alanine)']
    }

    # Use pandas for a clean, tabular output
    df = pd.DataFrame(mutation_data)
    print(df.to_string(index=False))

    print("\nSummary:")
    original_patch = "".join(mutation_data['Original AA (Code)'])
    mutated_patch = "AAAA" # All Alanine
    print(f"The site-directed mutagenesis experiment should change the sequence at positions 47-50")
    print(f"FROM: {mutation_data['Position'][0]}-{original_patch}-{mutation_data['Position'][-1]}")
    print(f"TO:   {mutation_data['Position'][0]}-{mutated_patch}-{mutation_data['Position'][-1]}")

# Execute the function to print the experimental design
design_mutagenesis_experiment()
<<<AAAA>>>