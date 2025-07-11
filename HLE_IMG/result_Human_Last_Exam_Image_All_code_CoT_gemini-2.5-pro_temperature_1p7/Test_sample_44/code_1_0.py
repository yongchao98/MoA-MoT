def solve_nmr_puzzle():
    """
    Analyzes 1H NMR data to identify the corresponding compound from a list of choices.
    The analysis is based on matching proton counts (total, aromatic) and aliphatic integration patterns.
    """
    # Step 1: Analyze the provided 1H NMR data
    nmr_data_raw = {
        8.19: {'integration': 1, 'type': 'aromatic'}, 7.79: {'integration': 1, 'type': 'aromatic'},
        7.47: {'integration': 1, 'type': 'aromatic'}, 7.38: {'integration': 1, 'type': 'aromatic'},
        6.98: {'integration': 1, 'type': 'aromatic'}, 6.63: {'integration': 1, 'type': 'aromatic'},
        6.61: {'integration': 1, 'type': 'aromatic'}, 4.19: {'integration': 4, 'type': 'aliphatic'},
        3.63: {'integration': 4, 'type': 'aliphatic'}, 3.21: {'integration': 2, 'type': 'aliphatic'},
        2.83: {'integration': 2, 'type': 'aliphatic'}, 1.98: {'integration': 2, 'type': 'aliphatic'}
    }

    # Calculate properties from NMR data
    nmr_aromatic_h = sum(v['integration'] for v in nmr_data_raw.values() if v['type'] == 'aromatic')
    nmr_aliphatic_integrations = sorted([v['integration'] for v in nmr_data_raw.values() if v['type'] == 'aliphatic'])
    nmr_total_h = sum(v['integration'] for v in nmr_data_raw.values())

    print("--- Step 1: Analysis of Provided 1H NMR Data ---")
    print(f"Total number of protons observed in the spectrum: {nmr_total_h}")
    print(f"Number of aromatic protons (from signals > 6.5 ppm): {nmr_aromatic_h}")
    print(f"Integration pattern of aliphatic protons (sorted): {nmr_aliphatic_integrations}\n")

    # Step 2: Define theoretical data for each compound
    # Protons from different parts of the base ligand structure:
    # Quinoline aromatic part: 3H
    # Tetrahydroquinoline aliphatic part (-CH2-CH2-CH2-): 6H (from 3 signals integrating to 2H each)
    # Piperazine part (-N-(CH2)2-N-): 8H (from 2 signals integrating to 4H each)
    # Thioamide NH proton: 1H (This is an acidic proton and is often not observed in the spectrum due to exchange)
    #
    # Substituents on the piperazine ring:
    # Pyridyl group: 4 aromatic H
    # Phenyl group: 5 aromatic H
    
    # Predictions will be made assuming the NH proton is NOT observed.
    
    compounds = {
        'A': {
            'description': 'Free ligand with a pyridyl group',
            'aromatic_H': 3 + 4,  # quinoline (3H) + pyridyl (4H)
            'aliphatic_integrations': sorted([2, 2, 2, 4, 4]) # 3xCH2 from tetrahydroquinoline + 2 sets of CH2 from piperazine
        },
        'C': {
            'description': 'Free ligand with a phenyl group',
            'aromatic_H': 3 + 5,  # quinoline (3H) + phenyl (5H)
            'aliphatic_integrations': sorted([2, 2, 2, 4, 4])
        },
        'B': {
            'description': 'cis-Zinc complex with two pyridyl ligands ([Zn(L_py)2])',
            'aromatic_H': 2 * (3 + 4),
            'aliphatic_integrations': sorted([2, 2, 2, 4, 4] * 2)
        },
        'D': {
            'description': 'trans-Zinc complex with two phenyl ligands ([Zn(L_ph)2])',
            'aromatic_H': 2 * (3 + 5),
            'aliphatic_integrations': sorted([2, 2, 2, 4, 4] * 2)
        },
        'E': {
            'description': 'trans-Zinc complex with two pyridyl ligands ([Zn(L_py)2])',
            'aromatic_H': 2 * (3 + 4),
            'aliphatic_integrations': sorted([2, 2, 2, 4, 4] * 2)
        }
    }
    
    for key, props in compounds.items():
        props['total_H'] = props['aromatic_H'] + sum(props['aliphatic_integrations'])

    print("--- Step 2 & 3: Comparing NMR Data to Each Structure ---")
    print("Predictions assume the single NH proton is not observed due to solvent exchange.\n")
    
    matching_compound = None
    for name, props in compounds.items():
        print(f"----- Checking Compound {name}: {props['description']} -----")
        print(f"  Predicted Total Protons: {props['total_H']}")
        print(f"  Predicted Aromatic Protons: {props['aromatic_H']}")
        print(f"  Predicted Aliphatic Integrations: {props['aliphatic_integrations']}")
        
        aromatic_match = (nmr_aromatic_h == props['aromatic_H'])
        aliphatic_match = (nmr_aliphatic_integrations == props['aliphatic_integrations'])
        total_match = (nmr_total_h == props['total_H'])

        if aromatic_match and aliphatic_match and total_match:
            print("  ==> Verdict: MATCH FOUND!")
            matching_compound = name
        else:
            print("  ==> Verdict: Does not match.")
        print("")

    # Step 4: Final conclusion
    print("--- Step 4: Final Conclusion ---")
    if matching_compound:
        print(f"The 1H NMR data perfectly matches the predicted spectrum for Compound {matching_compound}.")
        print("\nThe key matching points are:")
        print(f"1. Aromatic Proton Count: The NMR data shows {nmr_aromatic_h} aromatic protons. Compound {matching_compound} is predicted to have {compounds[matching_compound]['aromatic_H']} aromatic protons. This is a match.")
        print(f"2. Aliphatic Integration Pattern: The NMR data has aliphatic signals integrating to {nmr_aliphatic_integrations}. Compound {matching_compound} is predicted to have the same pattern. This is a match.")
        print(f"3. Total Proton Count: The total observed proton count is {nmr_total_h}, which matches the predicted count of {compounds[matching_compound]['total_H']} for Compound {matching_compound} (when the NH proton is excluded).")
        print("\nTherefore, the compound must be A.")
    else:
        print("No matching compound could be definitively identified based on the analysis.")

solve_nmr_puzzle()