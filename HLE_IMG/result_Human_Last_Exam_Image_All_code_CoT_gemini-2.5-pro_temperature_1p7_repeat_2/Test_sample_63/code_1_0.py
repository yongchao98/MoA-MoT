def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum to identify the correct molecular structure among candidates.
    This function prints the step-by-step reasoning.
    """
    
    print("### NMR Spectrum Analysis ###")
    print("The goal is to match the observed 1H NMR spectrum with one of the candidate structures (A-G, B-G, C-L, D-L).\n")

    # Observed data from the spectrum
    observed_signals = {
        "~1.2 ppm": "Triplet, 6H (corresponds to two equivalent -CH3 groups next to -CH2- groups)",
        "~2.3 ppm": "Singlet, 3H (corresponds to one -CH3 group with no proton neighbors)",
        "~2.8 ppm": "Quartet, 4H (corresponds to two equivalent -CH2- groups next to -CH3- groups)",
        "~3.6 ppm": "Singlet, 2H (corresponds to one -CH2- group with no proton neighbors)",
        "~7.1 ppm": "Multiplet, ~3H (corresponds to 3 protons on an aromatic ring)",
        "~9.0 ppm": "Broad Singlet, 1H (corresponds to an amide N-H proton)"
    }
    
    print("Step 1: Signals identified from the spectrum:")
    for shift, description in observed_signals.items():
        print(f"- Signal at {shift}: {description}")

    print("\nStep 2: Evaluating the candidates:")
    print(" - Candidates A-G and B-G can be eliminated. They lack the two ethyl groups (-CH2-CH3) that would produce the characteristic triplet (6H) and quartet (4H) seen in the spectrum.")
    print(" - Candidate D-L can be eliminated. It has two methyl groups on the aromatic ring, which would create a singlet of integration 6H at ~2.3 ppm. The spectrum shows a singlet of 3H.")
    
    print("\nStep 3: Verifying the best match (C-L):")
    
    structure_c_l_features = {
        "Amide N-H": "1H, broad singlet -> Matches the signal at ~9.0 ppm.",
        "Aromatic Protons": "3H, multiplet -> Matches the signal at ~7.1 ppm.",
        "Isolated Methylene (-CO-CH2-N)": "2H, singlet -> Matches the signal at ~3.6 ppm.",
        "Diethylamino group (-N(CH2CH3)2)": "Provides a 4H quartet (~2.8 ppm) and a 6H triplet (~1.2 ppm) -> Matches perfectly.",
        "Aromatic Methyl (-CH3)": "3H, singlet -> Matches the signal at ~2.3 ppm."
    }
    
    print("Candidate C-L has the following features:")
    for feature, explanation in structure_c_l_features.items():
        print(f"- {feature}: {explanation}")

    print("\nConclusion: The predicted spectrum for structure C-L is an excellent match for the observed spectrum.")

analyze_nmr_spectrum()