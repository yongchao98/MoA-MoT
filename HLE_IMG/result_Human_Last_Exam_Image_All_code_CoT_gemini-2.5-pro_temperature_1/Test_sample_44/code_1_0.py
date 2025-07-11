def solve_nmr_puzzle():
    """
    Analyzes 1H NMR data to identify the correct compound from a given list.
    """
    # 1. NMR Data Analysis
    nmr_integrations = [1, 1, 1, 1, 1, 1, 1, 4, 4, 2, 2, 2]
    nmr_shifts = [8.19, 7.79, 7.47, 7.38, 6.98, 6.63, 6.61, 4.19, 3.63, 3.21, 2.83, 1.98]
    total_protons_nmr = sum(nmr_integrations)

    aromatic_protons_nmr = sum(1 for shift in nmr_shifts if shift > 6.5)
    aliphatic_protons_nmr = total_protons_nmr - aromatic_protons_nmr

    print("--- NMR Data Analysis ---")
    print(f"Provided NMR integration values: {nmr_integrations}")
    print(f"Total number of protons from NMR data = {' + '.join(map(str, nmr_integrations))} = {total_protons_nmr}")
    print(f"Number of aromatic protons (δ > 6.5 ppm): {aromatic_protons_nmr}")
    print(f"Number of aliphatic protons (δ < 6.5 ppm): {aliphatic_protons_nmr}\n")

    # 2. Proton Count for Each Structure
    protons = {
        'tetrahydroquinoline': 9,  # 3 aromatic, 6 aliphatic
        'piperazine': 8,         # 8 aliphatic
        'pyridyl': 4,            # 4 aromatic
        'phenyl': 5,             # 5 aromatic
        'NH': 1
    }

    compounds = {
        'A': {
            'desc': '(Quinoline) + (NH) + (Piperazine) + (Pyridyl)',
            'total_H': protons['tetrahydroquinoline'] + protons['NH'] + protons['piperazine'] + protons['pyridyl'],
            'aromatic_H': 3 + 4, 'aliphatic_H': 6 + 8
        },
        'C': {
            'desc': '(Quinoline) + (NH) + (Piperazine) + (Phenyl)',
            'total_H': protons['tetrahydroquinoline'] + protons['NH'] + protons['piperazine'] + protons['phenyl'],
            'aromatic_H': 3 + 5, 'aliphatic_H': 6 + 8
        }
    }
    # Complexes have two ligands
    compounds['B'] = {'desc': 'Zn complex of ligand A', 'total_H': 2 * compounds['A']['total_H']}
    compounds['D'] = {'desc': 'Zn complex of ligand C', 'total_H': 2 * compounds['C']['total_H']}
    # Ligand for E has same proton count as A
    compounds['E'] = {'desc': 'Zn complex of a ligand similar to A', 'total_H': 2 * compounds['A']['total_H']}


    print("--- Structure Proton Count Analysis ---")
    print("Calculating the theoretical number of protons for each compound:")
    for name, data in compounds.items():
        print(f"Compound {name}: {data['desc']} = {data['total_H']} protons")

    # 3. Comparison and Hypothesis
    print("\n--- Comparison and Conclusion ---")
    print(f"The NMR data shows {total_protons_nmr} protons, which does not directly match any calculated total.")
    print("Hypothesis: The labile N-H proton is not observed in the spectrum.")
    print("Recalculating expected counts assuming the NH proton is not observed:")

    match_found = False
    for name in ['A', 'C']: # Only non-complexes are plausible
        data = compounds[name]
        adjusted_count = data['total_H'] - protons['NH']
        print(f"\nAdjusted count for Compound {name}: {data['total_H']} - 1 (NH) = {adjusted_count} protons.")
        
        if adjusted_count == total_protons_nmr:
            print(f"This matches the NMR total of {total_protons_nmr}.")
            print("Further validation by comparing proton types:")
            print(f"  - Compound {name} aromatic protons: {data['aromatic_H']}")
            print(f"  - NMR aromatic protons: {aromatic_protons_nmr}")
            print(f"  - Compound {name} aliphatic protons: {data['aliphatic_H']}")
            print(f"  - NMR aliphatic protons: {aliphatic_protons_nmr}")

            if data['aromatic_H'] == aromatic_protons_nmr and data['aliphatic_H'] == aliphatic_protons_nmr:
                print("\nThe distribution of aromatic and aliphatic protons also matches perfectly.")
                print(f"\nFinal Conclusion: The spectrum belongs to Compound {name}.")
                # The question asks for the Answer Choice letter. Compound A is option E.
                final_answer = 'E'
                match_found = True
                break
        else:
            print(f"This does not match the NMR total of {total_protons_nmr}.")

    if not match_found:
        print("\nCould not find a match.")
        
solve_nmr_puzzle()