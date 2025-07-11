def solve_nmr_puzzle():
    """
    Analyzes a 1H NMR spectrum to identify the correct molecular structure among candidates.
    This script codifies the logical steps of spectrum interpretation.
    """
    
    print("Analyzing the 1H NMR spectrum to find the best matching structure...\n")

    # Step 1: Identify key features from the spectrum.
    # The integration values are relative, normalized to the most downfield alkyl signal (triplet) being 6H.
    spectrum_features = {
        "diethyl_group_triplet": {"ppm": 1.2, "integration": 6, "multiplicity": "triplet"},
        "diethyl_group_quartet": {"ppm": 2.8, "integration": 4, "multiplicity": "quartet"},
        "amide_NH": {"ppm": 9.0, "integration": 1, "multiplicity": "broad singlet"},
        "linker_CH2": {"ppm": 3.4, "integration": 2, "multiplicity": "singlet"},
        "aryl_CH3": {"ppm": 2.3, "integration": 3, "multiplicity": "singlet"},
        "aromatic_H": {"ppm": 7.1, "integration": 3, "multiplicity": "complex multiplet"}
    }
    
    print("--- Step 1: Analysis of the Ethyl Group ---")
    t = spectrum_features["diethyl_group_triplet"]
    q = spectrum_features["diethyl_group_quartet"]
    print(f"The spectrum shows a triplet at ~{t['ppm']} ppm (integration = {t['integration']}H) and a quartet at ~{q['ppm']} ppm (integration = {q['integration']}H).")
    print("This pattern is the signature of a diethyl group [-N(CH2CH3)2].")
    print("This eliminates candidates A-G and B-G, which contain dimethyl groups.")
    print("Remaining candidates: C-L, D-L.\n")
    
    print("--- Step 2: Analysis of the Aryl Methyl Group ---")
    aryl_ch3_signal = spectrum_features["aryl_CH3"]
    print(f"The spectrum displays a singlet at ~{aryl_ch3_signal['ppm']} ppm with an integration of {aryl_ch3_signal['integration']}H.")
    print("This signal corresponds to one methyl group attached to the aromatic ring.\n")
    
    print("--- Step 3: Comparison with Remaining Candidates ---")
    print("Candidate C-L has one aryl methyl group, which would produce a singlet with 3H integration. This matches the spectrum.")
    print("Candidate D-L has two aryl methyl groups, which would produce a singlet with 6H integration. This contradicts the spectrum.")
    print("All other signals also match candidate C-L:")
    amide_s = spectrum_features["amide_NH"]
    linker_s = spectrum_features["linker_CH2"]
    print(f"  - Broad singlet at ~{amide_s['ppm']} ppm ({amide_s['integration']}H) -> Amide -NH")
    print(f"  - Singlet at ~{linker_s['ppm']} ppm ({linker_s['integration']}H) -> -CO-CH2-N-")
    print(f"  - Aromatic signals (~7.1 ppm, 3H) -> Trisubstituted benzene ring\n")
    
    print("--- Conclusion ---")
    print("The evidence overwhelmingly points to structure C-L as the correct one.")

solve_nmr_puzzle()
<<<C>>>