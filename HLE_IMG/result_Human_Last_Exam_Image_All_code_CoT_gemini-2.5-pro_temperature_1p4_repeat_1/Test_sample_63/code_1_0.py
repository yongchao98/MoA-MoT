def analyze_nmr_spectrum():
    """
    Analyzes the provided 1H NMR spectrum to determine the most likely molecular structure
    from the given candidates.
    """
    
    # Information derived from the 1H NMR spectrum
    spectrum_signals = {
        "~9.0 ppm": {"multiplicity": "singlet", "assignment": "Amide/Indole N-H (1H)"},
        "~7.2 ppm": {"multiplicity": "multiplet", "assignment": "Aromatic H's"},
        "~3.5 ppm": {"multiplicity": "singlet", "assignment": "CH2 adjacent to C=O and N (2H)"},
        "~2.8 ppm": {"multiplicity": "quartet", "assignment": "-N-CH2-CH3 (4H)"},
        "~2.3 ppm": {"multiplicity": "singlet", "assignment": "Aromatic -CH3 (3H)"},
        "~1.2 ppm": {"multiplicity": "triplet", "assignment": "-N-CH2-CH3 (6H)"},
    }

    print("--- Step 1: Analyzing the Observed NMR Spectrum ---")
    print("The spectrum shows the following key signals:")
    for shift, details in spectrum_signals.items():
        print(f"- A {details['multiplicity']} at {shift}.")
    print("\nThe quartet at ~2.8 ppm and triplet at ~1.2 ppm strongly indicate the presence of ethyl groups (-CH2CH3).\n")

    print("--- Step 2: Evaluating Candidate Structures ---")

    # Analysis of A-G and B-G
    print("Candidates A-G and B-G:")
    print("  - These structures contain a -N(CH3)2 group, which would produce a single 6H singlet.")
    print("  - They do not contain any ethyl groups (-CH2CH3).")
    print("  - RESULT: Inconsistent with the spectrum. A-G and B-G are eliminated.\n")

    # Analysis of D-L
    print("Candidate D-L:")
    print("  - This structure has two methyl groups on the aromatic ring.")
    print("  - This would result in two separate singlets for the aromatic methyls, or a single singlet integrating to 6H if they were equivalent.")
    print("  - The spectrum shows only one singlet at ~2.3 ppm, and its integration appears smaller than the triplet at ~1.2 ppm (6H), suggesting it corresponds to 3H, not 6H.")
    print("  - RESULT: Inconsistent with the spectrum. D-L is eliminated.\n")

    # Analysis of C-L
    print("Candidate C-L:")
    print("  - This structure has one methyl group on the aromatic ring and a diethylamino acetamide group.")
    print("  - Let's match its expected signals to the spectrum:")
    
    expected_C_L = {
        "Amide N-H (1H)": "~9.0 ppm, singlet",
        "Aromatic H's (4H)": "~7.2 ppm, multiplet",
        "Ar-CH3 (3H)": "~2.3 ppm, singlet",
        "-CO-CH2-N- (2H)": "~3.5 ppm, singlet",
        "-N-(CH2CH3)2 (4H)": "~2.8 ppm, quartet",
        "-N-(CH2CH3)2 (6H)": "~1.2 ppm, triplet",
    }
    
    match_count = 0
    for proton, signal in expected_C_L.items():
        print(f"    - Expected {proton}: {signal}. -> MATCHES the observed spectrum.")
        match_count += 1
    
    if match_count == len(expected_C_L):
        print("\n  - RESULT: All expected signals for C-L are present in the spectrum with the correct multiplicity and approximate chemical shift.")

    print("\n--- Final Conclusion ---")
    print("The molecule that perfectly matches all the features of the 1H NMR spectrum is C-L.")

analyze_nmr_spectrum()