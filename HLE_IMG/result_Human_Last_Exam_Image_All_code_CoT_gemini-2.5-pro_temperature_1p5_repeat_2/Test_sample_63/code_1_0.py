def analyze_nmr_spectrum():
    """
    Analyzes an 1H NMR spectrum and compares it to four candidate drug structures
    to find the best match.
    """
    # Step 1: Define the data extracted from the experimental NMR spectrum
    spectrum_data = {
        '9.0 ppm': {'protons': 1, 'type': 'broad singlet', 'assignment': 'Amide N-H'},
        '7.0-7.5 ppm': {'protons': '4 (approx.)', 'type': 'multiplet', 'assignment': 'Aromatic protons'},
        '3.5 ppm': {'protons': 2, 'type': 'singlet', 'assignment': '-CO-CH2-N-'},
        '2.8 ppm': {'protons': 4, 'type': 'quartet', 'assignment': '-N-(CH2CH3)2'},
        '2.4 ppm': {'protons': 3, 'type': 'singlet', 'assignment': 'Aryl-CH3'},
        '1.2 ppm': {'protons': 6, 'type': 'triplet', 'assignment': '-N-(CH2CH3)2'}
    }

    print("--- Analysis of the Experimental 1H NMR Spectrum ---")
    for shift, data in spectrum_data.items():
        print(f"Signal at ~{shift}: {data['protons']}H, {data['type']}. Corresponds to: {data['assignment']}.")
    print("-" * 50)

    # Step 2: Define the predicted signals for each candidate structure
    predictions = {
        'A-G': "Indole N-H (1H, s), Aromatic H (5H, m), -CH2-CH2- (4H, m), -N(CH3)2 (6H, s). Mismatch: No ethyl group (quartet/triplet) signals.",
        'B-G': "Indole N-H (1H, s), Aromatic H (5H, m), -CH2- (2H, s), -N(CH3)2 (6H, s). Mismatch: No ethyl group signals; expects a 6H singlet.",
        'C-L': "Amide N-H (1H, s), Aromatic H (4H, m), Aryl-CH3 (3H, s), -CO-CH2-N- (2H, s), -N(CH2CH3)2 giving a 4H quartet and a 6H triplet. Perfect Match.",
        'D-L': "Amide N-H (1H, s), Aromatic H (3H, m), Aryl-CH3 (6H, s from two methyls), -CO-CH2-N- (2H, s), -N(CH2CH3)2 giving a 4H quartet and a 6H triplet. Mismatch: Expects a 6H singlet for Aryl-CH3, but spectrum shows a 3H singlet."
    }

    print("--- Evaluating Candidate Structures ---")
    for molecule, prediction in predictions.items():
        print(f"Candidate {molecule}:")
        print(f"  Predicted Signals & Comparison: {prediction}\n")

    print("--- Conclusion ---")
    print("Structure C-L is the only candidate whose predicted 1H NMR signals for chemical shift, integration (number of protons), and splitting pattern (multiplicity) all align perfectly with the provided spectrum.")

analyze_nmr_spectrum()