def analyze_nmr_spectrum():
    """
    Analyzes the provided 1H NMR spectrum to identify the correct molecular structure.
    """
    print("Step 1: Initial analysis of the NMR spectrum for key patterns.")
    print("The spectrum shows a distinct quartet (q) around 2.8 ppm and a triplet (t) around 1.2 ppm.")
    print("This q-t pattern is a classic signature of an ethyl group (-CH2-CH3).")
    print("-" * 30)

    print("Step 2: Eliminate candidates that do not have an ethyl group.")
    print("Candidates A-G and B-G have N,N-dimethyl groups, not ethyl groups.")
    print("Therefore, candidates A-G and B-G are eliminated.")
    print("The possible structures are C-L and D-L, as both contain N,N-diethyl groups.")
    print("-" * 30)

    print("Step 3: Differentiate between C-L and D-L based on integration.")
    print("The key difference is the number of methyl groups on the aromatic ring:")
    print(" - Structure C-L has ONE aromatic methyl group (Ar-CH3).")
    print(" - Structure D-L has TWO aromatic methyl groups (Ar-(CH3)2).")
    print("\nThis will affect the integration of the aromatic methyl singlet (~2.3 ppm) relative to the ethyl methyl triplet (~1.2 ppm).")
    print("-" * 30)
    
    print("Step 4: Predict and compare the integration ratios.")
    # Predicted proton counts for the key signals
    protons_c_l = {'Aromatic Methyls (singlet at ~2.3 ppm)': 3, 'Ethyl Methyls (triplet at ~1.2 ppm)': 6}
    protons_d_l = {'Aromatic Methyls (singlet at ~2.3 ppm)': 6, 'Ethyl Methyls (triplet at ~1.2 ppm)': 6}

    print("Prediction for C-L:")
    print(f" - Protons for Ar-CH3 (singlet): {protons_c_l['Aromatic Methyls (singlet at ~2.3 ppm)']}H")
    print(f" - Protons for ethyl -CH3 (triplet): {protons_c_l['Ethyl Methyls (triplet at ~1.2 ppm)']}H")
    ratio_c_l = protons_c_l['Aromatic Methyls (singlet at ~2.3 ppm)'] / protons_c_l['Ethyl Methyls (triplet at ~1.2 ppm)']
    print(f" - Expected integration ratio (Ar-CH3 singlet / Ethyl-CH3 triplet) = {protons_c_l['Aromatic Methyls (singlet at ~2.3 ppm)']}:{protons_c_l['Ethyl Methyls (triplet at ~1.2 ppm)']} or {ratio_c_l:.2f}")
    print()

    print("Prediction for D-L:")
    print(f" - Protons for Ar-(CH3)2 (singlet): {protons_d_l['Aromatic Methyls (singlet at ~2.3 ppm)']}H")
    print(f" - Protons for ethyl -CH3 (triplet): {protons_d_l['Ethyl Methyls (triplet at ~1.2 ppm)']}H")
    ratio_d_l = protons_d_l['Aromatic Methyls (singlet at ~2.3 ppm)'] / protons_d_l['Ethyl Methyls (triplet at ~1.2 ppm)']
    print(f" - Expected integration ratio (Ar-CH3 singlet / Ethyl-CH3 triplet) = {protons_d_l['Aromatic Methyls (singlet at ~2.3 ppm)']}:{protons_d_l['Ethyl Methyls (triplet at ~1.2 ppm)']} or {ratio_d_l:.2f}")
    print("-" * 30)
    
    print("Step 5: Conclusion from spectral evidence.")
    print("By visual inspection of the spectrum, the integrated area of the triplet at ~1.2 ppm is clearly larger than the area of the singlet at ~2.3 ppm.")
    print("The area of the triplet appears to be about double the area of the singlet, matching the 1:2 ratio (3H:6H) predicted for C-L.")
    print("The 1:1 ratio predicted for D-L is inconsistent with the spectrum.")
    print("\nTherefore, the most possible structure is C-L.")

analyze_nmr_spectrum()