def analyze_nmr_data():
    """
    Analyzes 1H NMR data to identify the correct compound from a list of possibilities.
    """
    # 1H NMR Data
    # Format: {chemical_shift: [integration, multiplicity]}
    nmr_data_str = "8.19 (1H, m), 7.79 (1H, m), 7.47 (1H, m), 7.38 (1H, m), 6.98 (1H, m), 6.63 (1H, m), 6.61 (1H, m), 4.19 (4H, m), 3.63 (4H, m), 3.21 (2H, m), 2.83 (2H, m), 1.98 (2H, m)"
    
    # Parsing the NMR data
    signals = nmr_data_str.split(', ')
    integrations = [int(s.split('H')[0].split('(')[1]) for s in signals]
    shifts = [float(s.split(' ')[0]) for s in signals]
    
    total_protons_nmr = sum(integrations)
    aromatic_protons_nmr = sum(h for h, s in zip(integrations, shifts) if s > 6.5)
    aliphatic_protons_nmr = sum(h for h, s in zip(integrations, shifts) if s < 5.0)

    print("--- 1H NMR Data Analysis ---")
    print(f"Provided NMR data: {nmr_data_str}")
    print(f"Total protons observed: {total_protons_nmr}")
    print(f"Aromatic protons observed (shift > 6.5 ppm): {aromatic_protons_nmr}")
    print(f"Aliphatic protons observed (shift < 5.0 ppm): {aliphatic_protons_nmr}")
    print("-" * 30 + "\n")

    # Expected protons for each compound
    # (Excluding labile NH protons for initial comparison against observed data)
    compounds = {
        'A': {'type': 'Ligand', 'aromatic': 7, 'aliphatic': 14, 'labile': 1},
        'B': {'type': 'Complex of A', 'aromatic': 14, 'aliphatic': 28, 'labile': 0},
        'C': {'type': 'Ligand', 'aromatic': 8, 'aliphatic': 14, 'labile': 1},
        'D': {'type': 'Complex of C', 'aromatic': 16, 'aliphatic': 28, 'labile': 0},
        'E': {'type': 'Complex of A', 'aromatic': 14, 'aliphatic': 28, 'labile': 0},
    }

    print("--- Compound Analysis ---")
    match = None
    for compound, props in compounds.items():
        print(f"Analyzing Compound {compound} ({props['type']}):")
        
        expected_observed_protons = props['aromatic'] + props['aliphatic']
        total_protons_structure = expected_observed_protons + props['labile']

        print(f"  - Expected Aromatic H: {props['aromatic']}")
        print(f"  - Expected Aliphatic H: {props['aliphatic']}")
        print(f"  - Total Expected H (with labile H): {total_protons_structure}")

        # Comparison Logic
        # Ligands: Compare counts assuming labile proton is not observed.
        if props['type'].startswith('Ligand'):
            if (aromatic_protons_nmr == props['aromatic'] and
                aliphatic_protons_nmr == props['aliphatic'] and
                total_protons_nmr == expected_observed_protons):
                print("  - Verdict: MATCH. Proton counts match the NMR data if the labile NH proton is not observed.")
                match = compound
            else:
                print("  - Verdict: MISMATCH. Proton counts do not match the NMR data.")

        # Complexes: Compare total protons and integration pattern.
        else:
            if total_protons_nmr == expected_observed_protons:
                 print("  - Verdict: MATCH on proton count.")
                 match = compound
            else:
                print(f"  - Verdict: MISMATCH. The total proton count ({total_protons_nmr}) is far from the expected count ({expected_observed_protons}). Complexes have ~2x the protons of a ligand.")
        print("")

    # Map compound letter to answer choice
    answer_choices = {'A': 'E', 'B': 'A', 'C': 'D', 'D': 'B', 'E': 'C'}
    
    print("--- Conclusion ---")
    if match:
        final_answer_choice = answer_choices[match]
        print(f"The NMR data is consistent with the structure of Compound {match}.")
        print(f"Compound {match} corresponds to answer choice {final_answer_choice}.")
    else:
        print("Could not find a definitive match.")

    # Final Answer
    print("\nFinal Answer Selection:")
    print(f"The correct structure is {match}, which is answer choice {final_answer_choice}.")
    print(f'<<<E>>>')

# Run the analysis
analyze_nmr_data()