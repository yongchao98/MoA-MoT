import re

def solve_nmr_puzzle():
    """
    Analyzes 1H NMR data to identify the corresponding compound from a list of choices.
    """
    # 1. Analyze the 1H NMR Data
    nmr_data_str = "8.19 (1H, m), 7.79 (1H, m), 7.47 (1H, m), 7.38 (1H, m), 6.98 (1H, m), 6.63 (1H, m), 6.61 (1H, m), 4.19 (4H, m), 3.63 (4H, m), 3.21 (2H, m), 2.83 (2H, m), 1.98 (2H, m)"

    # Parse the string to get shift and integration
    signals = re.findall(r'([\d\.]+) \((\d+)H, m\)', nmr_data_str)
    
    nmr_total_h = 0
    nmr_aromatic_h = 0
    nmr_aliphatic_h = 0

    print("Analyzing 1H NMR Data:")
    print("Signal(ppm)\tIntegration(#H)\tRegion")
    print("------------------------------------------")
    for shift_str, integration_str in signals:
        shift = float(shift_str)
        integration = int(integration_str)
        nmr_total_h += integration
        
        if shift > 6.0:
            region = "Aromatic"
            nmr_aromatic_h += integration
        else:
            region = "Aliphatic"
            nmr_aliphatic_h += integration
        print(f"{shift:.2f}\t\t{integration}\t\t{region}")
    
    print("\n--- NMR Data Summary ---")
    print(f"Total Protons (H) from NMR: {nmr_aromatic_h} + {nmr_aliphatic_h} = {nmr_total_h}")
    print(f"Aromatic Protons (H): {nmr_aromatic_h}")
    print(f"Aliphatic Protons (H): {nmr_aliphatic_h}")
    print("--------------------------\n")

    # 2. Analyze the Candidate Compounds
    # Counts are based on the structures.
    # For A and C (ligands), the single NH proton is assumed to not be observed.
    # For B, D, E (complexes), there are two ligands per complex.
    compounds = {
        'A': {
            # 3 (quinoline arom) + 4 (pyridine arom) = 7 aromatic H
            # 6 (quinoline aliph) + 8 (piperazine aliph) = 14 aliphatic H
            'aromatic': 7, 'aliphatic': 14, 'total': 21,
            'description': "Ligand with pyridyl group"
        },
        'B': {
            # 2 ligands * (7 aromatic H + 14 aliphatic H)
            'aromatic': 14, 'aliphatic': 28, 'total': 42,
            'description': "Zinc complex with two 'A' ligands"
        },
        'C': {
            # 3 (quinoline arom) + 5 (phenyl arom) = 8 aromatic H
            # 6 (quinoline aliph) + 8 (piperazine aliph) = 14 aliphatic H
            'aromatic': 8, 'aliphatic': 14, 'total': 22,
            'description': "Ligand with phenyl group"
        },
        'D': {
            # 2 ligands * (8 aromatic H + 14 aliphatic H)
            'aromatic': 16, 'aliphatic': 28, 'total': 44,
            'description': "Zinc complex with two 'C' ligands"
        },
        'E': {
            # Same as B, just a different coordination isomer drawing
            # 2 ligands * (7 aromatic H + 14 aliphatic H)
            'aromatic': 14, 'aliphatic': 28, 'total': 42,
            'description': "Zinc complex isomer of B"
        }
    }
    
    # 3. Compare and Conclude
    print("Comparing NMR Data with Candidate Compounds:")
    match_found = None
    for name, props in compounds.items():
        print(f"\n--- Checking Compound {name} ({props['description']}) ---")
        
        compound_total_h = props['total']
        compound_aromatic_h = props['aromatic']
        compound_aliphatic_h = props['aliphatic']
        
        # We present the comparison as an "equation" of counts
        print(f"NMR counts:         Aromatic={nmr_aromatic_h}, Aliphatic={nmr_aliphatic_h}, Total={nmr_total_h}")
        print(f"Compound {name} counts: Aromatic={compound_aromatic_h}, Aliphatic={compound_aliphatic_h}, Total={compound_total_h}")

        if (nmr_total_h == compound_total_h and
            nmr_aromatic_h == compound_aromatic_h and
            nmr_aliphatic_h == compound_aliphatic_h):
            print("\nResult: PERFECT MATCH!")
            match_found = name
        else:
            print("\nResult: No Match.")

    # 4. Final Answer Mapping
    # The question provides a mapping from answer choice to compound label
    answer_choices = {'A': 'B', 'B': 'D', 'C': 'E', 'D': 'C', 'E': 'A'}

    final_answer_choice = None
    for choice, compound_label in answer_choices.items():
        if compound_label == match_found:
            final_answer_choice = choice
            break

    if final_answer_choice:
        print("\n==============================================")
        print(f"Conclusion: The NMR data corresponds to Compound {match_found}.")
        print(f"The correct answer choice for Compound {match_found} is {final_answer_choice}.")
        print("==============================================")
        # Final answer format required by the user prompt
        print(f"\n<<<{final_answer_choice}>>>")
    else:
        print("\nCould not find a matching compound.")

# Run the analysis
solve_nmr_puzzle()