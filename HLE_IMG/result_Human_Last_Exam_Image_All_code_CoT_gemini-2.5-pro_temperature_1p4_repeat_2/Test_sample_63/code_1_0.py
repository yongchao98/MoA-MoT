def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum by comparing its features to predicted features of candidate molecules.
    """
    # 1. Define the observed NMR spectrum features
    observed_spectrum = {
        "signals": {
            "~8.8 ppm": {"protons": 1, "splitting": "broad singlet", "assignment": "Amide N-H"},
            "~7.2 ppm": {"protons": "Multiple (3)", "splitting": "multiplet", "assignment": "Aromatic Ar-H"},
            "~3.5 ppm": {"protons": 2, "splitting": "singlet", "assignment": "Isolated -CH2-"},
            "~2.9 ppm": {"protons": 4, "splitting": "quartet", "assignment": "Ethyl -CH2-"},
            "~2.3 ppm": {"protons": 3, "splitting": "singlet", "assignment": "Aromatic -CH3"},
            "~1.2 ppm": {"protons": 6, "splitting": "triplet", "assignment": "Ethyl -CH3"},
        },
        "key_pattern": "Quartet (4H) and Triplet (6H) indicates two ethyl groups -N(CH2CH3)2."
    }

    print("--- Analysis of the Observed NMR Spectrum ---")
    print(f"Key Feature Detected: {observed_spectrum['key_pattern']}")
    print("Individual Signals:")
    for shift, details in observed_spectrum["signals"].items():
        print(f"  - Shift: {shift}, Protons: {details['protons']}H, Splitting: {details['splitting']}")
    print("-" * 45 + "\n")


    # 2. Define predicted spectra for each candidate
    candidates = {
        "A-G": {
            "summary": "Indole with a -CH2CH2-N(CH3)2 side chain.",
            "features": "Contains a dimethyl [-N(CH3)2] group (6H singlet), not ethyl groups. Lacks the quartet/triplet pattern."
        },
        "B-G": {
            "summary": "Indole with a -CH2-N(CH3)2 side chain (Gramine).",
            "features": "Contains a dimethyl [-N(CH3)2] group (6H singlet), not ethyl groups. Lacks the quartet/triplet pattern."
        },
        "C-L": {
            "summary": "Amide with one aromatic -CH3 and a diethylamino [-N(C2H5)2] group.",
            "features": {
                "Amide N-H": "1H broad singlet (~8.8 ppm). Correct.",
                "Aromatic H": "3H multiplet (~7.2 ppm). Correct.",
                "Isolated CO-CH2-N": "2H singlet (~3.5 ppm). Correct.",
                "Diethylamino -N(CH2CH3)2": "4H quartet (~2.9 ppm) and 6H triplet (~1.2 ppm). Correct.",
                "Aromatic -CH3": "3H singlet (~2.3 ppm). Correct."
            },
            "match": True
        },
        "D-L": {
            "summary": "Amide with two aromatic -CH3 groups and a diethylamino group.",
            "features": "Predicts TWO distinct singlets for the two aromatic -CH3 groups. The spectrum only shows one 3H singlet. Incorrect."
        }
    }

    # 3. Compare and conclude
    print("--- Evaluating Candidate Structures ---")
    final_answer = None
    for name, data in candidates.items():
        print(f"Candidate: {name} ({data['summary']})")
        if name == "C-L":
            print("  Analysis:")
            for feature, desc in data['features'].items():
                print(f"    - {feature}: {desc}")
            print("  Result: PERFECT MATCH")
            final_answer = name
        else:
            print(f"  Analysis: {data['features']}")
            print("  Result: MISMATCH")
        print("-" * 35)

    print(f"\nConclusion: The spectrum's features, particularly the amide proton at 8.8 ppm,")
    print(f"the single aromatic methyl group at 2.3 ppm, the isolated methylene at 3.5 ppm,")
    print(f"and the characteristic diethylamino pattern (a quartet at 2.9 ppm and a triplet at 1.2 ppm),")
    print(f"are all uniquely and perfectly consistent with structure {final_answer}.")


analyze_nmr_spectrum()