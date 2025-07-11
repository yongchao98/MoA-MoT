def analyze_nmr_spectrum():
    """
    This function analyzes the given 1H NMR spectrum and matches it to the
    most plausible chemical structure among the candidates.
    """

    print("Step-by-step Analysis of the 1H NMR Spectrum for Structure C-L:")
    print("================================================================")
    
    # Define the observed signals and their assignments to structure C-L
    signals = {
        "Amide N-H": "~9.0 ppm (broad singlet, 1H)",
        "Aromatic Protons (4H)": "~7.0-7.5 ppm (multiplet, 4H)",
        "Methylene -CO-CH2-N-": "~3.4 ppm (singlet, 2H)",
        "Methylene in Ethyl groups -N-(CH2CH3)2": "~2.8 ppm (quartet, 4H)",
        "Aromatic Methyl -CH3": "~2.3 ppm (singlet, 3H)",
        "Methyl in Ethyl groups -N-(CH2CH3)2": "~1.1 ppm (triplet, 6H)",
    }
    
    # Print the rationale for each signal based on structure C-L
    for group, spec in signals.items():
        print(f"- Signal for {group}: {spec}")

    print("\nJustification for choice:")
    print("1. The presence of a 6H triplet (~1.1 ppm) and a 4H quartet (~2.8 ppm) is a clear indication of two equivalent ethyl groups (-C2H5). This matches the diethylamino group in C-L and D-L, and rules out A-G and B-G which have dimethylamino groups.")
    print("2. The spectrum shows a single 3H singlet at ~2.3 ppm, which corresponds to the one methyl group on the aromatic ring of structure C-L.")
    print("3. Structure D-L is ruled out because it has two methyl groups on its aromatic ring, which would lead to signals different from the single 3H singlet observed.")
    
    print("\nFinal Conclusion:")
    print("The spectrum unambiguously matches structure C-L.")

analyze_nmr_spectrum()