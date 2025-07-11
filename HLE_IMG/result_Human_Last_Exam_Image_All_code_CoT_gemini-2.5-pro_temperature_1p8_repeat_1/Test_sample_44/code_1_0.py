import pandas as pd

def analyze_nmr_data():
    """
    This function analyzes the provided 1H NMR data to determine the structure.
    """
    # Step 1: Define and analyze the 1H NMR data
    nmr_data = {
        'shift': [8.19, 7.79, 7.47, 7.38, 6.98, 6.63, 6.61, 4.19, 3.63, 3.21, 2.83, 1.98],
        'integration': [1, 1, 1, 1, 1, 1, 1, 4, 4, 2, 2, 2]
    }
    
    total_protons_nmr = sum(nmr_data['integration'])
    aromatic_protons_nmr = sum(h for s, h in zip(nmr_data['shift'], nmr_data['integration']) if s > 6.0)
    aliphatic_protons_nmr = sum(h for s, h in zip(nmr_data['shift'], nmr_data['integration']) if s < 6.0)

    print("--- 1H NMR Data Analysis ---")
    print(f"Total protons observed in NMR spectrum: {total_protons_nmr}")
    print(f"Aromatic protons (shift > 6.0 ppm) in NMR spectrum: {aromatic_protons_nmr}")
    print(f"Aliphatic protons (shift < 6.0 ppm) in NMR spectrum: {aliphatic_protons_nmr}")
    print("-" * 30 + "\n")

    # Step 2: Define proton counts for each compound
    # L-py is the ligand in A, B, E. L-ph is the ligand in C, D.
    compounds = {
        'A': {
            'description': 'Ligand with pyridine substituent',
            'aromatic_H': 7,  # 3 (quinoline) + 4 (pyridine)
            'aliphatic_H': 14, # 6 (tetrahydroquinoline) + 8 (piperazine)
            'exchangeable_H': 1 # NH
        },
        'C': {
            'description': 'Ligand with phenyl substituent',
            'aromatic_H': 8,  # 3 (quinoline) + 5 (phenyl)
            'aliphatic_H': 14, # 6 (tetrahydroquinoline) + 8 (piperazine)
            'exchangeable_H': 1 # NH
        },
        'B': {
            'description': 'Zn complex with two ligands of type A',
            'aromatic_H': 2 * 7,
            'aliphatic_H': 2 * 14,
            'exchangeable_H': 0
        },
        'D': {
            'description': 'Zn complex with two ligands of type C',
            'aromatic_H': 2 * 8,
            'aliphatic_H': 2 * 14,
            'exchangeable_H': 0
        },
        'E': {
            'description': 'Zn complex with two ligands of type A (identical to B)',
            'aromatic_H': 2 * 7,
            'aliphatic_H': 2 * 14,
            'exchangeable_H': 0
        }
    }
    
    # Pre-calculate total protons for each compound
    for name, props in compounds.items():
        props['total_H'] = props['aromatic_H'] + props['aliphatic_H'] + props['exchangeable_H']
        props['total_H_no_exchange'] = props['aromatic_H'] + props['aliphatic_H']

    # Step 3: Compare NMR data with each compound
    print("--- Comparison with Compounds ---")
    
    match_found = False
    for name, props in compounds.items():
        print(f"\nAnalyzing Compound {name}: {props['description']}")
        
        print(f"  - Expected total protons: {props['total_H']}")
        print(f"  - Expected total protons (if NH not observed): {props['total_H_no_exchange']}")
        print(f"  - Expected aromatic protons: {props['aromatic_H']}")
        print(f"  - Expected aliphatic protons: {props['aliphatic_H']}")
        
        # Check for match, assuming NH proton might not be observed
        if (props['total_H_no_exchange'] == total_protons_nmr and
            props['aromatic_H'] == aromatic_protons_nmr and
            props['aliphatic_H'] == aliphatic_protons_nmr):
            
            print("\n>>> MATCH FOUND <<<")
            print(f"The NMR data is consistent with Compound {name}, assuming the exchangeable NH proton is not observed.")
            print(f"Comparison Summary for Compound {name}:")
            print(f"  - Total Protons:  Observed={total_protons_nmr}, Expected(w/o NH)={props['total_H_no_exchange']}")
            print(f"  - Aromatic Protons: Observed={aromatic_protons_nmr}, Expected={props['aromatic_H']}")
            print(f"  - Aliphatic Protons: Observed={aliphatic_protons_nmr}, Expected={props['aliphatic_H']}")
            match_found = True
            final_answer_compound = name
            break # Stop after finding the first match
        else:
            print("  - Result: No match.")
    
    if not match_found:
        print("\nCould not find a matching compound.")
    else:
        print("\n--- Final Conclusion ---")
        print(f"The 1H NMR data corresponds to compound {final_answer_compound}.")
        print("This is because the total number of protons (21), the number of aromatic protons (7), and the number of aliphatic protons (14) perfectly match the structure of compound A, assuming the single NH proton is not observed due to solvent exchange, a common phenomenon in NMR spectroscopy.")
        # Find the multiple choice letter for the compound
        choices = {'A': 'E', 'B': 'A', 'C': 'D', 'D': 'B', 'E': 'C'}
        print(f"Compound {final_answer_compound} corresponds to multiple choice option {choices[final_answer_compound]}.")


if __name__ == '__main__':
    analyze_nmr_data()