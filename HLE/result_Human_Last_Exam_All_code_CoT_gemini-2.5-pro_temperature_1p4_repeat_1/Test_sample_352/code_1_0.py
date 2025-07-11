import pandas as pd

def design_mutagenesis_experiment():
    """
    Outlines the plan and justification for a site-directed mutagenesis experiment
    to relieve an autoinhibitory negatively charged patch in a protein.
    """

    # Define the original and proposed mutant amino acids at the specified positions
    original_patch = {
        "Position": [47, 48, 49, 50],
        "Original Amino Acid": ["S (Serine)", "E (Glutamate)", "E (Glutamate)", "D (Aspartate)"],
        "Reason for Change": [
            "Phosphorylation site, becomes highly negative",
            "Negatively charged (acidic)",
            "Negatively charged (acidic)",
            "Negatively charged (acidic)"
        ]
    }

    mutant_patch = {
        "Position": [47, 48, 49, 50],
        "Proposed Mutant": ["A (Alanine)", "A (Alanine)", "A (Alanine)", "A (Alanine)"],
        "Justification": [
            "Cannot be phosphorylated, removes negative charge potential",
            "Neutral, removes negative charge",
            "Neutral, removes negative charge",
            "Neutral, removes negative charge"
        ]
    }

    # Create dataframes for better visualization
    original_df = pd.DataFrame(original_patch).set_index("Position")
    mutant_df = pd.DataFrame(mutant_patch).set_index("Position")

    print("Experiment: Relieve Autoinhibitory Negative Charge in Protein X")
    print("=" * 65)
    print("\nThe original patch at positions 47-50 is S-E-E-D and is highly negatively charged:")
    print(original_df)

    print("\nProposed Solution: Replace the entire patch with Alanine (A-A-A-A)")
    print("-" * 65)
    print("This comprehensively neutralizes the patch. Here is the justification for each mutation:")
    print(mutant_df)

    print("\nFinal Proposed Mutant Sequence at positions 47-50: A-A-A-A")
    print("This corresponds to the mutations: S47A, E48A, E49A, D50A.")
    print("This is the best replacement to test the hypothesis that the negative")
    print("charge of the patch is responsible for autoinhibition.")

# Execute the function to print the plan
design_mutagenesis_experiment()