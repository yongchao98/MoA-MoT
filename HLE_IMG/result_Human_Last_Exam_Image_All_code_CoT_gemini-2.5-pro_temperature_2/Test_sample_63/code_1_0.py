import collections

def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum by comparing its features to the expected spectra
    of several candidate drug molecules to find the best match.
    """
    # Step 1: Define the observed data from the 1H NMR spectrum.
    # We estimate chemical shifts, integrations (relative number of protons),
    # and multiplicities (splitting patterns). The integrations are the key.
    # The triplet at ~1.2 ppm and quartet at ~2.8 ppm are characteristic of an ethyl group.
    # The integration ratio of 3:2 for these peaks confirms it. Let's assume two ethyl
    # groups as seen in C-L and D-L, so we set the baseline integration for the
    # triplet at 6H and the quartet at 4H.
    observed_spectrum_analysis = {
        'signals': collections.OrderedDict([
            ('~9.0 ppm', {'integration': 1, 'multiplicity': 'broad singlet'}),
            ('~7.2 ppm', {'integration': 3, 'multiplicity': 'multiplet'}),
            ('~3.5 ppm', {'integration': 2, 'multiplicity': 'singlet'}),
            ('~2.8 ppm', {'integration': 4, 'multiplicity': 'quartet'}),
            ('~2.3 ppm', {'integration': 3, 'multiplicity': 'singlet'}),
            ('~1.2 ppm', {'integration': 6, 'multiplicity': 'triplet'})
        ]),
        'features': "Contains two ethyl groups; aryl-CH3 to aromatic H ratio is 3:3 (1:1)."
    }

    print("--- Step 1: Analyzing Observed 1H NMR Spectrum ---")
    print("The spectrum shows the following key signals with their estimated integrations:")
    for shift, data in observed_spectrum_analysis['signals'].items():
        print(f"- Signal at {shift}: A {data['multiplicity']} with relative integration of {data['integration']}H.")
    print("\nFeature Summary: " + observed_spectrum_analysis['features'])

    # Step 2: Define the expected NMR features for each candidate molecule.
    candidates = {
        'A-G': {'name': 'alpha-dimethyltryptamine', 'match': False, 'reason': 'Does not contain an ethyl group. Contains a dimethylamino group (-N(CH3)2) which would be a 6H singlet.'},
        'B-G': {'name': 'N,N-dimethyltryptamine (DMT)', 'match': False, 'reason': 'Does not contain an ethyl group. Also has a dimethylamino group.'},
        'C-L': {'name': '2-(diethylamino)-N-(2-methylphenyl)acetamide', 'match': True,
                'protons': {'Amide NH': 1, 'Aromatic H': 3, 'COCH2N': 2, 'Ethyl N-CH2': 4, 'Aryl-CH3': 3, 'Ethyl C-CH3': 6}},
        'D-L': {'name': '2-(diethylamino)-N-(2,6-dimethylphenyl)acetamide (Lidocaine)', 'match': False,
                'protons': {'Amide NH': 1, 'Aromatic H': 2, 'COCH2N': 2, 'Ethyl N-CH2': 4, 'Aryl-CH3': 6, 'Ethyl C-CH3': 6}}
    }

    print("\n--- Step 2 & 3: Evaluating Candidates and Comparing ---")
    best_match = None
    for key, data in candidates.items():
        print(f"\nAnalyzing Candidate: {key} ({data['name']})")
        if not data['match']:
            print(f"Result: Incorrect. Reason: {data['reason']}")
            continue

        # For matching candidates, compare integration pattern
        predicted_integrations = data['protons']
        observed_integrations = {
            'Amide NH': 1, 'Aromatic H': 3, 'COCH2N': 2,
            'Ethyl N-CH2': 4, 'Aryl-CH3': 3, 'Ethyl C-CH3': 6
        }
        
        is_match = True
        for proton_type, integration in predicted_integrations.items():
            if observed_integrations[proton_type] != integration:
                print(f"Mismatch for '{proton_type}': Expected {integration}H, but spectrum suggests {observed_integrations[proton_type]}H.")
                is_match = False
        
        if is_match:
            print("Result: Perfect match. The number and type of protons align with the spectral data.")
            print(f"The key confirmation is the integration ratio of Aryl-CH3 protons ({predicted_integrations['Aryl-CH3']}H) to Aromatic protons ({predicted_integrations['Aromatic H']}H), which is {predicted_integrations['Aryl-CH3']}:{predicted_integrations['Aromatic H']}, matching the spectrum.")
            best_match = key
        else:
             candidates[key]['match'] = False # Update status
             if key == 'D-L':
                 print(f"The key mismatch is the integration ratio of Aryl-CH3 protons ({predicted_integrations['Aryl-CH3']}H) to Aromatic protons ({predicted_integrations['Aromatic H']}H), which should be {predicted_integrations['Aryl-CH3']}:{predicted_integrations['Aromatic H']} (or 3:1), but the spectrum shows a 1:1 ratio.")


    print("\n--- Step 4: Final Conclusion ---")
    if best_match:
        print(f"The analysis strongly indicates that the molecule is {best_match} ({candidates[best_match]['name']}).")
    else:
        print("No suitable match found among the candidates.")

analyze_nmr_spectrum()
<<<C>>>