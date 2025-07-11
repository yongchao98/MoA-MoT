def analyze_nmr_spectrum():
    """
    Analyzes an 1H NMR spectrum to identify the correct molecular structure
    from a list of candidates.
    """
    print("--- Step 1: Analysis of the 1H NMR Spectrum ---")
    print("The spectrum shows the following signals:")
    print("Signal 1: Chemical Shift = ~9.0 ppm, Multiplicity = broad singlet, Integration = 1H. Assignment: Amide N-H proton.")
    print("Signal 2: Chemical Shift = ~7.1 ppm, Multiplicity = multiplet. Assignment: Aromatic protons (Ar-H).")
    print("Signal 3: Chemical Shift = ~3.4 ppm, Multiplicity = singlet, Integration = 2H. Assignment: Methylene group (-CH2-) adjacent to carbonyl and nitrogen.")
    print("Signal 4: Chemical Shift = ~2.8 ppm, Multiplicity = quartet, Integration = 4H. Assignment: Two equivalent methylene groups (-CH2-) next to methyl groups.")
    print("Signal 5: Chemical Shift = ~2.3 ppm, Multiplicity = singlet, Integration = 3H. Assignment: Aromatic methyl group (Ar-CH3).")
    print("Signal 6: Chemical Shift = ~1.1 ppm, Multiplicity = triplet, Integration = 6H. Assignment: Two equivalent methyl groups (-CH3) next to methylene groups.")
    print("\nKey Feature: The quartet (4H) and triplet (6H) strongly indicate the presence of two equivalent ethyl groups: -N(CH2CH3)2.\n")

    print("--- Step 2: Evaluation of Candidate Structures ---")

    # Candidate A-G and B-G
    print("\nAnalysis of A-G and B-G:")
    print("Expected: Both structures have a dimethylamino group, -N(CH3)2.")
    print("This would produce a 6H singlet.")
    print("Observed: The spectrum has a 6H triplet and a 4H quartet, not a 6H singlet.")
    print("Conclusion: A-G and B-G are incorrect.\n")

    # Candidate D-L
    print("Analysis of D-L:")
    print("Expected: This structure has two methyl groups on the aromatic ring.")
    print("This would produce a signal for aromatic methyls integrating to 6H.")
    print("Observed: The spectrum shows a singlet at ~2.3 ppm integrating to only 3H.")
    print("Conclusion: D-L is incorrect.\n")

    # Candidate C-L
    print("Analysis of C-L:")
    print("This structure is 2-(diethylamino)-N-(2-methylphenyl)acetamide.")
    print("Let's check its expected signals against the observed spectrum:")
    print(f"- Amide N-H (1H): Matches the broad singlet at ~9.0 ppm.")
    print(f"- Aromatic Protons (3H): Matches the multiplet at ~7.1 ppm.")
    print(f"- Aromatic Methyl (3H): Matches the singlet at ~2.3 ppm.")
    print(f"- -C(=O)-CH2-N- (2H): This methylene has no neighbors, predicting a singlet. Matches the 2H singlet at ~3.4 ppm.")
    print(f"- Diethylamino -N(CH2CH3)2: Predicts a 4H quartet for the -CH2- and a 6H triplet for the -CH3. This perfectly matches the quartet at ~2.8 ppm and the triplet at ~1.1 ppm.")
    print("Conclusion: All signals in the spectrum are consistent with structure C-L.\n")

    print("--- Final Conclusion ---")
    print("Based on the detailed analysis, the NMR spectrum corresponds to structure C-L.")

analyze_nmr_spectrum()
<<<C>>>