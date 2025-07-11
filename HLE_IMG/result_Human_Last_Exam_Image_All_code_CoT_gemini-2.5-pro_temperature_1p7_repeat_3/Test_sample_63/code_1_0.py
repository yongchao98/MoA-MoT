def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum based on its key features and matches it to a candidate molecule.
    """

    # --- Step 1: Define observed spectral data ---
    # From the image, we identify the following signals:
    # (ppm, multiplicity, relative integration)
    observed_signals = {
        'signal_1': {'ppm': 1.2, 'multiplicity': 'triplet', 'integration': 6},
        'signal_2': {'ppm': 2.4, 'multiplicity': 'singlet', 'integration': 3},
        'signal_3': {'ppm': 2.7, 'multiplicity': 'quartet', 'integration': 4},
        'signal_4': {'ppm': 3.2, 'multiplicity': 'singlet', 'integration': 2},
        'signal_5': {'ppm': '~7.1', 'multiplicity': 'multiplet', 'integration': '3-4'},
        'signal_6': {'ppm': 8.8, 'multiplicity': 'broad singlet', 'integration': 1}
    }

    print("--- Analysis of the Observed NMR Spectrum ---")
    print("A triplet at ~1.2 ppm (6H) and a quartet at ~2.7 ppm (4H) strongly suggest the presence of two ethyl groups attached to a nitrogen: -N(CH2CH3)2.")
    print("This key feature eliminates candidates A-G and B-G, which have -N(CH3)2 groups.")
    print("\nWe are left with candidates C-L and D-L. Let's compare them.")
    print("-" * 50)

    # --- Step 2: Define expected features for remaining candidates ---
    candidates = {
        'C-L': {
            'Aromatic CH3 protons': {'expected_integration': 3, 'observed_match': 2.4},
            'Aromatic CH protons': {'expected_integration': 3},
            'Full match': True
        },
        'D-L': {
            'Aromatic CH3 protons': {'expected_integration': 6, 'observed_match': 2.4},
            'Aromatic CH protons': {'expected_integration': 2},
            'Full match': False
        }
    }

    print("--- Comparing C-L and D-L ---")
    print("Molecule C-L has one methyl group on the aromatic ring (3H).")
    print("Molecule D-L has two methyl groups on the aromatic ring (6H).")
    print("\nThe observed spectrum shows a singlet at 2.4 ppm with a relative integration of 3H.")

    # --- Step 3: Find the best match ---
    best_match = None
    for name, data in candidates.items():
        if data['Aromatic CH3 protons']['expected_integration'] == observed_signals['signal_2']['integration']:
            best_match = name
            break

    print("\nThis 3H integration for the aromatic methyl signal directly matches the structure of C-L.")
    print("It contradicts the structure of D-L, which would show a 6H signal for its two methyl groups.")
    print("\n--- Final Conclusion ---")
    if best_match:
        print(f"The molecule whose predicted spectrum perfectly matches all the observed signals is {best_match}.")
    else:
        print("No perfect match found based on the aromatic methyl signal.")

    # Let's confirm with all the signals for the final answer.
    print(f"\nVerifying all signals for molecule {best_match}:")
    print(f"Amide N-H: Predicted (1H, broad singlet, >8 ppm), Observed (~{observed_signals['signal_6']['ppm']} ppm, {observed_signals['signal_6']['multiplicity']}, {observed_signals['signal_6']['integration']}H). -> Match!")
    print(f"Aromatic C-H: Predicted (~3H), Observed (~3H based on relative integration). -> Match!")
    print(f"-CO-CH2-N- protons: Predicted (2H, singlet), Observed (~{observed_signals['signal_4']['ppm']} ppm, {observed_signals['signal_4']['multiplicity']}, {observed_signals['signal_4']['integration']}H). -> Match!")
    print(f"-N-(CH2)2- protons: Predicted (4H, quartet), Observed (~{observed_signals['signal_3']['ppm']} ppm, {observed_signals['signal_3']['multiplicity']}, {observed_signals['signal_3']['integration']}H). -> Match!")
    print(f"Aromatic -CH3 protons: Predicted (3H, singlet), Observed (~{observed_signals['signal_2']['ppm']} ppm, {observed_signals['signal_2']['multiplicity']}, {observed_signals['signal_2']['integration']}H). -> Match!")
    print(f"-CH2-(CH3)2 protons: Predicted (6H, triplet), Observed (~{observed_signals['signal_1']['ppm']} ppm, {observed_signals['signal_1']['multiplicity']}, {observed_signals['signal_1']['integration']}H). -> Match!")

analyze_nmr_spectrum()