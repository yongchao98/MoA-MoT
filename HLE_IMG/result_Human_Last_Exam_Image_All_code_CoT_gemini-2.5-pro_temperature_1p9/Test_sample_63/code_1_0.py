def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum by comparing its features to predicted spectra of candidate molecules.
    """
    # Step 1: Define the signals observed in the experimental NMR spectrum.
    # Each signal is a dictionary with its properties.
    experimental_spectrum = [
        {'ppm': 8.8, 'integration': 1, 'splitting': 'broad singlet', 'description': 'likely an N-H proton'},
        {'ppm': 7.2, 'integration': 3, 'splitting': 'multiplet', 'description': 'protons on a substituted benzene ring'},
        {'ppm': 3.4, 'integration': 2, 'splitting': 'singlet', 'description': 'a CH2 group with no adjacent protons'},
        {'ppm': 2.7, 'integration': 4, 'splitting': 'quartet', 'description': 'a CH2 group next to a CH3 group'},
        {'ppm': 2.3, 'integration': 3, 'splitting': 'singlet', 'description': 'a CH3 group with no adjacent protons'},
        {'ppm': 1.1, 'integration': 6, 'splitting': 'triplet', 'description': 'two CH3 groups next to a CH2 group'},
    ]

    print("--- Analysis of the Experimental NMR Spectrum ---")
    total_protons = 0
    for signal in experimental_spectrum:
        print(f"Signal at ~{signal['ppm']} ppm: Integration={signal['integration']}H, Splitting={signal['splitting']}. Notes: {signal['description']}.")
        total_protons += signal['integration']
    print(f"\nTotal proton integration from spectrum: {total_protons}H")

    # Step 2: Identify key structural patterns from the spectrum.
    print("\n--- Identifying Key Structural Fragments ---")
    has_ethyl_group = False
    quartet_4h = any(s['integration'] == 4 and s['splitting'] == 'quartet' for s in experimental_spectrum)
    triplet_6h = any(s['integration'] == 6 and s['splitting'] == 'triplet' for s in experimental_spectrum)

    if quartet_4h and triplet_6h:
        has_ethyl_group = True
        print("1. Found a 4H quartet and a 6H triplet. This is a classic signature for a diethylamino group: -N(CH2CH3)2.")
    else:
        print("1. No clear pattern for an ethyl group was found.")

    aromatic_methyl_3h = any(s['integration'] == 3 and s['splitting'] == 'singlet' and 2.0 < s['ppm'] < 2.5 for s in experimental_spectrum)
    if aromatic_methyl_3h:
        print("2. Found a 3H singlet around 2.3 ppm, characteristic of a methyl group on a benzene ring (Ar-CH3).")

    # Step 3: Define and evaluate candidate structures.
    candidates = {
        "A-G": {
            "protons": 16,
            "has_diethyl_group": False,
            "aromatic_methyls_H": 0
        },
        "B-G": {
            "protons": 16,
            "has_diethyl_group": False,
            "aromatic_methyls_H": 0
        },
        "C-L": {
            "protons": 19,
            "has_diethyl_group": True,
            "aromatic_methyls_H": 3,
            "description": "Has one aromatic methyl group (3H), three aromatic protons (3H), one amide N-H (1H), one -CO-CH2- group (2H), and one -N(CH2CH3)2 group (10H)."
        },
        "D-L": {
            "protons": 22,
            "has_diethyl_group": True,
            "aromatic_methyls_H": 6,
            "description": "Has two aromatic methyl groups (6H total)."
        }
    }

    print("\n--- Comparing with Candidate Structures ---")
    best_match = "None"
    for name, properties in candidates.items():
        match_score = 0
        print(f"\nEvaluating Candidate: {name}")

        # Check for diethyl group
        if properties["has_diethyl_group"] == has_ethyl_group:
            print(f"- Match: Presence of diethyl group is consistent with the spectrum.")
            match_score += 1
        else:
            print(f"- Mismatch: Presence of diethyl group is NOT consistent with the spectrum.")
            continue # In this case, it's a critical mismatch

        # Check total proton count
        if properties["protons"] == total_protons:
             print(f"- Match: Total proton count ({properties['protons']}H) matches the spectrum's integration ({total_protons}H).")
             match_score += 1
        else:
             print(f"- Mismatch: Total proton count ({properties['protons']}H) does NOT match spectrum ({total_protons}H).")

        # Check aromatic methyl protons
        if aromatic_methyl_3h and properties["aromatic_methyls_H"] == 3:
            print(f"- Match: The 3H aromatic methyl signal is consistent with this structure.")
            match_score += 1
        elif not aromatic_methyl_3h and properties["aromatic_methyls_H"] == 0:
            match_score += 1 # Consistent if neither has it
        else:
            print(f"- Mismatch: The spectrum shows a 3H singlet for an aromatic methyl, but this structure predicts {properties['aromatic_methyls_H']}H.")


        if match_score == 3: # A perfect match on our key criteria
            best_match = name
            print(f"\nConclusion for {name}: This structure is a perfect match for all key features of the spectrum.")
            print(f"Detailed breakdown for {best_match}: {properties['description']}")


    # Final Conclusion
    print("\n--- FINAL CONCLUSION ---")
    if best_match != "None":
        print(f"The molecule that best fits the NMR data is {best_match}.")
    else:
        print("Could not find a perfect match among the candidates.")

analyze_nmr_spectrum()
<<<C>>>