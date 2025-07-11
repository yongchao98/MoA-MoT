import pandas as pd

def analyze_lec_stability():
    """
    Analyzes the expected stability of LECs based on three different Ir(III) complexes.
    """

    # Data representing the features of each complex relevant to device stability
    data = {
        'Complex': ['Complex 1', 'Complex 2', 'Complex 3'],
        'Cyclometalating Ligand': ['ppy (standard)', 'ppy (standard)', 'dfppy (fluorinated)'],
        'Ancillary Ligand': ['bpy (standard)', 'Large, planar ligand', 'dtbbpy (bulky groups)'],
        'Expected Aggregation': ['Moderate', 'High', 'Low'],
        'Chemical Stability Enhancement': ['None', 'None', 'Fluorination (High)'],
        'Morphological Stability Enhancement': ['None', 'None', 'Bulky Groups (High)']
    }

    # Assign scores based on features to quantify stability
    # Stability Score = (Chemical Stability) + (Morphological Stability) - (Aggregation Tendency)
    scores = {
        'None': 0, 'Moderate': 0, 'High': -1, 'Low': 1, 'Fluorination (High)': 1, 'Bulky Groups (High)': 1
    }

    df = pd.DataFrame(data)

    df['Stability Score'] = (df['Chemical Stability Enhancement'].map(scores) +
                             df['Morphological Stability Enhancement'].map(scores) +
                             df['Expected Aggregation'].map(scores))

    print("--- Analysis of Ir(III) Complexes for LEC Stability ---\n")
    print("The operational stability of a Light-Emitting Electrochemical Cell (LEC) is heavily")
    print("influenced by the molecular design of the emissive complex. Key factors include:")
    print("  - Chemical Stability: Resistance to degradation under electrical stress.")
    print("  - Morphological Stability: Resistance to aggregation and crystallization in the thin film.\n")

    print("Let's evaluate each complex:\n")

    # Print analysis for each complex
    for index, row in df.iterrows():
        complex_name = row['Complex']
        score = row['Stability Score']
        print(f"--- {complex_name} ---")
        if complex_name == 'Complex 1':
            print("This is a baseline reference complex. It has moderate stability but is prone to gradual degradation.")
        elif complex_name == 'Complex 2':
            print("Features a large, planar ligand that promotes strong intermolecular aggregation (π-π stacking).")
            print("This leads to luminescence quenching and poor morphological stability, reducing device lifetime.")
        elif complex_name == 'Complex 3':
            print("This complex is designed for high stability:")
            print("  1. Fluorination: The C-F bonds enhance chemical and electrochemical stability.")
            print("  2. Bulky tert-butyl groups: These prevent aggregation, improving solid-state emission and film morphology.")
        print(f"Calculated Stability Score: {score}\n")

    # Determine the most stable complex
    most_stable_idx = df['Stability Score'].idxmax()
    most_stable_complex = df.loc[most_stable_idx]

    print("--- Conclusion ---")
    print(f"Based on the analysis, {most_stable_complex['Complex']} achieves the highest stability score.")
    print(f"Its design directly addresses common degradation pathways by incorporating both")
    print(f"fluorination for chemical stability and bulky groups for morphological stability.")
    print("\nTherefore, LECs based on Complex 3 are expected to be the most stable.")

analyze_lec_stability()
<<<C>>>