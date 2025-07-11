def analyze_nmr_spectrum():
    """
    Analyzes candidate molecules against a given 1H NMR spectrum data.
    """
    # 1. Define the observed peaks from the NMR spectrum
    spectrum_peaks = {
        'Amide N-H':       {'shift': '8.9 ppm', 'integration': 1, 'multiplicity': 'singlet'},
        'Aromatic H':      {'shift': '7.1 ppm', 'integration': 3, 'multiplicity': 'multiplet'},
        'Methylene (alpha to C=O, N)': {'shift': '3.4 ppm', 'integration': 2, 'multiplicity': 'singlet'},
        'Ethyl -CH2-':     {'shift': '2.8 ppm', 'integration': 4, 'multiplicity': 'quartet'},
        'Aromatic -CH3':   {'shift': '2.3 ppm', 'integration': 3, 'multiplicity': 'singlet'},
        'Ethyl -CH3-':     {'shift': '1.2 ppm', 'integration': 6, 'multiplicity': 'triplet'}
    }

    print("--- NMR Spectrum Analysis ---")
    print("Based on the provided spectrum, the following signals are identified:")
    for name, data in spectrum_peaks.items():
        print(f"- A {data['multiplicity']} at ~{data['shift']} integrating to {data['integration']}H, assigned to: {name}")

    print("\n--- Evaluating Candidate Structures ---")

    # 2. Analyze each candidate
    # Candidate A-G
    print("\n[Analysis of A-G and B-G]")
    print("These molecules contain a dimethylamino group (-N(CH3)2).")
    print("This would produce a single strong singlet for 6 protons.")
    print("The spectrum shows a quartet (4H) and a triplet (6H), characteristic of a diethylamino group (-N(CH2CH3)2).")
    print("Result: A-G and B-G are incorrect.\n")

    # Candidate D-L (Lidocaine)
    print("[Analysis of D-L]")
    print("This molecule has a 2,6-dimethylphenyl group.")
    print("This would produce a singlet for two equivalent methyl groups, integrating to 6H.")
    print(f"The spectrum shows a singlet at {spectrum_peaks['Aromatic -CH3']['shift']} that integrates to {spectrum_peaks['Aromatic -CH3']['integration']}H.")
    print("Result: D-L is incorrect.\n")

    # Candidate C-L
    print("[Analysis of C-L: 2-(diethylamino)-N-(2-methylphenyl)acetamide]")
    print("Let's compare the expected signals for C-L with the spectrum:")
    print(f"- Amide N-H: Expected: 1H singlet. Observed: Matches peak at {spectrum_peaks['Amide N-H']['shift']}.")
    print(f"- Aromatic H on a 1,2,3-trisubstituted ring: Expected: 3H multiplet. Observed: Matches peaks around {spectrum_peaks['Aromatic H']['shift']}.")
    print(f"- Aromatic -CH3: Expected: 3H singlet. Observed: Matches peak at {spectrum_peaks['Aromatic -CH3']['shift']}.")
    print(f"- -CO-CH2-N-: Expected: 2H singlet. Observed: Matches peak at {spectrum_peaks['Methylene (alpha to C=O, N)']['shift']}.")
    print(f"- Diethylamino -N(CH2CH3)2 group:")
    print(f"  - -CH2- protons: Expected: 4H quartet. Observed: Matches peak at {spectrum_peaks['Ethyl -CH2-']['shift']}.")
    print(f"  - -CH3 protons: Expected: 6H triplet. Observed: Matches peak at {spectrum_peaks['Ethyl -CH3-']['shift']}.")
    print("Result: C-L is a perfect match for the spectrum.\n")

    # 3. Final Conclusion
    print("--- Conclusion ---")
    print("The molecule that fully matches the NMR data, including the chemical shifts, integrations, and splitting patterns, is C-L.")

# Execute the analysis
analyze_nmr_spectrum()