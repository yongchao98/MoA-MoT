import pandas as pd

def analyze_nmr():
    """
    This function analyzes the NMR spectrum and candidate structures to find the correct match.
    """
    print("Step 1: Analyzing the provided 1H NMR Spectrum")
    print("--------------------------------------------------")
    
    # Data extracted by observing the image of the NMR spectrum
    observed_signals = {
        '~8.8 ppm': {'Protons': '1H', 'Splitting': 'Broad Singlet', 'Assignment': 'Amide or N-H on a ring'},
        '~7.1 ppm': {'Protons': '~3H', 'Splitting': 'Multiplet', 'Assignment': 'Aromatic Protons (Ar-H)'},
        '~3.5 ppm': {'Protons': '2H', 'Splitting': 'Singlet', 'Assignment': 'Isolated -CH2- group'},
        '~2.8 ppm': {'Protons': '4H', 'Splitting': 'Quartet', 'Assignment': 'Ethyl group -CH2-'},
        '~2.3 ppm': {'Protons': '6H', 'Splitting': 'Singlet', 'Assignment': 'Two equivalent -CH3 groups'},
        '~1.2 ppm': {'Protons': '6H', 'Splitting': 'Triplet', 'Assignment': 'Ethyl group -CH3'},
    }
    
    # Creating a DataFrame for better visualization
    df_observed = pd.DataFrame(observed_signals).T
    print(df_observed)
    
    print("\nKey features from the spectrum:")
    print("- A pair of signals (triplet and quartet) characteristic of two equivalent ethyl groups [-N(CH2CH3)2].")
    print("- A broad signal far downfield (~8.8 ppm), typical for an amide N-H proton.")
    print("- A large singlet at ~2.3 ppm, suggesting two equivalent methyl groups (6H).")

    print("\nStep 2: Predicting NMR Spectra for Candidate Structures")
    print("--------------------------------------------------------")
    
    predictions = {
        'A-G': 'Has a -N(CH3)2 group (dimethylamino), not -N(C2H5)2 (diethylamino). No ethyl group signals expected. --> Incorrect.',
        'B-G': 'Has a -N(CH3)2 group, not -N(C2H5)2. No ethyl group signals expected. --> Incorrect.',
        'C-L': 'Represents a Lidocaine-like structure. It contains a diethylamino group [-N(C2H5)2], which would produce a triplet (~6H) and a quartet (~4H). It also has an amide N-H bond, which explains the broad singlet at ~8.8 ppm. The structure as drawn shows one Ar-CH3 group, but the spectrum\'s large singlet at 2.3 ppm strongly suggests two equivalent methyl groups (6H), which matches the actual structure of Lidocaine (N-(2,6-dimethylphenyl)acetamide). This is the best fit despite the drawing inaccuracy.',
        'D-L': 'This structure has a diethylamino group, but it also has a second, different ethyl group attached to the aromatic ring. This would result in two different sets of triplet/quartet signals. The spectrum only shows one set. --> Incorrect.'
    }
    
    print("\nPredictions for each candidate:")
    for candidate, analysis in predictions.items():
        print(f"- {candidate}: {analysis}")

    print("\nStep 3: Conclusion")
    print("------------------")
    print("Based on the analysis:")
    print("1. Candidates A-G and B-G are ruled out because they contain dimethylamino groups, while the spectrum clearly shows signals for a diethylamino group.")
    print("2. Candidate D-L is ruled out because it would show signals for two distinct ethyl groups, but the spectrum only shows one set.")
    print("3. Candidate C-L is the only structure containing the correct core fragments: a diethylamino group and an amide N-H. The signals at ~1.2 ppm (triplet, 6H), ~2.8 ppm (quartet, 4H), ~3.5 ppm (singlet, 2H for CO-CH2-N), and ~8.8 ppm (singlet, 1H for N-H) are all perfectly consistent with this structure type. The strong singlet at ~2.3 ppm indicates two equivalent methyl groups on the aromatic ring, making the structure Lidocaine.")
    
    print("\nTherefore, the spectrum corresponds to the structure class represented by C-L.")


# Run the analysis
analyze_nmr()
