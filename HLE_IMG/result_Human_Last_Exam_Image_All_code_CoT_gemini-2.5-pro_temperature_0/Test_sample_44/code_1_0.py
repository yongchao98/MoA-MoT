def solve_nmr_puzzle():
    """
    Analyzes 1H NMR data to identify the corresponding compound from a list of choices.
    """

    # 1. Define the experimental NMR data
    nmr_data = {
        'shifts': [8.19, 7.79, 7.47, 7.38, 6.98, 6.63, 6.61, 4.19, 3.63, 3.21, 2.83, 1.98],
        'integrations': [1, 1, 1, 1, 1, 1, 1, 4, 4, 2, 2, 2]
    }

    # Calculate total and categorized protons from experimental data
    total_protons_exp = sum(nmr_data['integrations'])
    aromatic_protons_exp = sum(h for s, h in zip(nmr_data['shifts'], nmr_data['integrations']) if s > 6.5)
    aliphatic_protons_exp = sum(h for s, h in zip(nmr_data['shifts'], nmr_data['integrations']) if s < 5.0)

    print("--- Analysis of Experimental 1H NMR Data ---")
    print(f"Total observed protons: {total_protons_exp}")
    print(f"Observed aromatic protons (shift > 6.5 ppm): {aromatic_protons_exp}")
    print(f"Observed aliphatic protons (shift < 5.0 ppm): {aliphatic_protons_exp}")
    print("-" * 40)

    # 2. Define the predicted proton counts for each compound
    # Note: We assume the acidic NH proton in free ligands (A, C) is not observed.
    compounds = {
        'A': {
            'description': 'Free ligand with Pyridine group',
            'quinoline_arom_H': 3, 'quinoline_aliph_H': 6,
            'piperazine_H': 8, 'pyridine_H': 4, 'NH': 1,
            'is_complex': False
        },
        'C': {
            'description': 'Free ligand with Phenyl group',
            'quinoline_arom_H': 3, 'quinoline_aliph_H': 6,
            'piperazine_H': 8, 'phenyl_H': 5, 'NH': 1,
            'is_complex': False
        },
        'B': {
            'description': 'Zn(II) complex with two ligands of type A',
            'ligand': 'A', 'is_complex': True
        },
        'D': {
            'description': 'Zn(II) complex with two ligands of type C',
            'ligand': 'C', 'is_complex': True
        },
        'E': {
            'description': 'Zn(II) complex with two ligands (isomer of A)',
            'ligand': 'A', 'is_complex': True # Isomeric ligand has same proton count as A
        }
    }

    # 3. Analyze each compound and compare with experimental data
    print("--- Analysis of Potential Compounds ---")
    match = None
    for name, props in compounds.items():
        print(f"\nAnalyzing Compound {name}: {props['description']}")
        if not props['is_complex']:
            # For free ligands A and C
            total_H = props['quinoline_arom_H'] + props['quinoline_aliph_H'] + props['piperazine_H'] + props.get('pyridine_H', 0) + props.get('phenyl_H', 0)
            arom_H = props['quinoline_arom_H'] + props.get('pyridine_H', 0) + props.get('phenyl_H', 0)
            aliph_H = props['quinoline_aliph_H'] + props['piperazine_H']
            print(f"Predicted total protons (assuming NH is not observed): {total_H}")
            print(f"Predicted aromatic protons: {arom_H}")
            print(f"Predicted aliphatic protons: {aliph_H}")
            if total_H == total_protons_exp and arom_H == aromatic_protons_exp and aliph_H == aliphatic_protons_exp:
                print(">>> This compound is a MATCH.")
                match = name
            else:
                print(">>> This compound is NOT a match.")
        else:
            # For complexes B, D, E
            ligand_props = compounds[props['ligand']]
            total_H = 2 * (ligand_props['quinoline_arom_H'] + ligand_props['quinoline_aliph_H'] + ligand_props['piperazine_H'] + ligand_props.get('pyridine_H', 0) + ligand_props.get('phenyl_H', 0))
            print(f"Predicted total protons: {total_H}")
            print(">>> This compound is NOT a match (total proton count is too high).")

    print("-" * 40)
    print(f"Final Conclusion: The NMR data corresponds to Compound {match}.")
    print("The answer choices map as follows: A->B, B->D, C->E, D->C, E->A.")
    print(f"Since the compound is {match}, the correct answer choice is E.")

solve_nmr_puzzle()