def analyze_nmr_spectrum():
    """
    Analyzes an 1H NMR spectrum to determine the most likely chemical structure among candidates.
    This function will print the step-by-step reasoning.
    """

    print("Step 1: Analyze the key signals in the 1H NMR spectrum.")
    spectrum_signals = {
        "~8.9 ppm": "Broad Singlet (1H)",
        "~7.1 ppm": "Multiplet (Aromatic Region)",
        "~3.2 ppm": "Singlet",
        "~2.6 ppm": "Quartet",
        "~2.3 ppm": "Singlet",
        "~1.1 ppm": "Triplet"
    }
    for shift, description in spectrum_signals.items():
        print(f"- Signal at {shift}: {description}")

    print("\nStep 2: Identify structural fragments from the spectrum.")
    print("The combination of a Triplet at ~1.1 ppm and a Quartet at ~2.6 ppm is a classic pattern for an ethyl group (-CH2CH3).")
    print("This strongly suggests the molecule contains one or more ethyl groups.")

    print("\nStep 3: Evaluate the candidate structures.")
    print("- Structures A-G and B-G do not contain ethyl groups. They contain dimethylamino groups (-N(CH3)2), which would appear as a single 6H singlet. These structures are inconsistent with the spectrum.")
    print("- Structures C-L and D-L both contain two ethyl groups attached to a nitrogen atom [-N(C2H5)2], which would produce a 6H triplet and a 4H quartet. These are the plausible candidates.")

    print("\nStep 4: Differentiate between Structure C-L and D-L.")
    print("The key difference is the number of protons on the aromatic ring.")
    print("- Structure C-L (1,2-disubstituted ring) has 4 aromatic protons.")
    print("- Structure D-L (1,2,4-trisubstituted ring) has 3 aromatic protons.")
    print("By comparing the integration of the aromatic multiplet (~7.1 ppm) with other signals (e.g., the 4H quartet at ~2.6 ppm), the area corresponds to 4 protons.")
    print("This evidence strongly favors Structure C-L.")

    print("\nStep 5: Final assignment of signals for Structure C-L.")
    final_assignment = {
        "Structure C-L Proton": "Predicted Signal (Integration)",
        "Amide H (-NH-)": "~8.9 ppm (broad singlet, 1H)",
        "Aromatic H (Ar-H)": "~7.1 ppm (multiplet, 4H)",
        "Methylene H (-CO-CH2-)": "~3.2 ppm (singlet, 2H)",
        "Ethyl Methylene H (-N-(CH2)2)": "~2.6 ppm (quartet, 4H)",
        "Aromatic Methyl H (Ar-CH3)": "~2.3 ppm (singlet, 3H)",
        "Ethyl Methyl H (-(CH2CH3)2)": "~1.1 ppm (triplet, 6H)"
    }
    print("="*60)
    for proton, signal in final_assignment.items():
        print(f"{proton:<30} -> {signal}")
    print("="*60)
    print("\nConclusion: The spectrum is a perfect match for structure C-L.")

# Execute the analysis
analyze_nmr_spectrum()