def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum by comparing its features to predicted spectra of candidate molecules.
    """
    # Observed features from the 1H NMR spectrum
    observed_spectrum = {
        "~9.0 ppm": {"protons": 1, "type": "broad singlet", "group": "Amide N-H"},
        "~7.2 ppm": {"protons": 3, "type": "multiplet", "group": "Aromatic C-H"},
        "~3.4 ppm": {"protons": 2, "type": "singlet", "group": "-CO-CH2-N-"},
        "~2.8 ppm": {"protons": 4, "type": "quartet", "group": "Ethyl -CH2-"},
        "~2.5 ppm": {"protons": 3, "type": "singlet", "group": "Aryl -CH3"},
        "~1.2 ppm": {"protons": 6, "type": "triplet", "group": "Ethyl -CH3"}
    }
    
    # Predicted features for each candidate molecule
    candidates = {
        "A-G": "Contains N(CH3)2 (6H singlet), lacks ethyl group. -> Incorrect.",
        "B-G": "Contains N(CH3)2 (6H singlet), lacks ethyl group. -> Incorrect.",
        "C-L": {
            "Amide N-H": "1H, broad singlet -> Match",
            "Aromatic C-H": "3H, multiplet -> Match",
            "Aryl -CH3": "3H, singlet -> Match",
            "-CO-CH2-N-": "2H, singlet -> Match",
            "Ethyl group -N(CH2CH3)2": "4H quartet + 6H triplet -> Match"
        },
        "D-L": {
            "Amide N-H": "1H, broad singlet -> Match",
            "Aromatic C-H": "2H, singlet -> Mismatch (spectrum suggests ~3H)",
            "Aryl -CH3": "6H (2xCH3), singlet -> Mismatch (spectrum suggests ~3H)",
            "-CO-CH2-N-": "2H, singlet -> Match",
            "Ethyl group -N(CH2CH3)2": "4H quartet + 6H triplet -> Match"
        }
    }
    
    print("### NMR Spectrum Analysis ###\n")
    print("1. Key features identified from the spectrum:")
    print("- A quartet (~2.8 ppm) and a triplet (~1.2 ppm) indicate the presence of ethyl groups (-CH2CH3).")
    print("- This immediately rules out candidates A-G and B-G, which lack ethyl groups.\n")
    
    print("2. Comparing remaining candidates C-L and D-L:")
    
    print("\n--- Analysis for Candidate C-L ---")
    for feature, desc in candidates["C-L"].items():
        print(f"- {feature}: {desc}")
    print("Conclusion: All features of candidate C-L match the observed spectrum.")
    
    print("\n--- Analysis for Candidate D-L ---")
    for feature, desc in candidates["D-L"].items():
        print(f"- {feature}: {desc}")
    print("Conclusion: Candidate D-L mismatches on the number of aromatic protons (2H vs ~3H) and aryl methyl protons (6H vs ~3H).\n")
    
    print("### Final Decision ###")
    print("The spectrum is most consistent with structure C-L.")

analyze_nmr_spectrum()