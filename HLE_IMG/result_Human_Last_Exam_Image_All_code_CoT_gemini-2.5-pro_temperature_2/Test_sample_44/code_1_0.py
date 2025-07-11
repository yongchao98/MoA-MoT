def identify_compound_from_nmr():
    """
    Analyzes 1H NMR data to identify the corresponding chemical structure
    from a given set of options.
    """
    # 1. Define the provided 1H NMR data
    # Format: (chemical_shift_ppm, integration_H)
    nmr_data = [
        (8.19, 1), (7.79, 1), (7.47, 1), (7.38, 1), (6.98, 1), (6.63, 1), (6.61, 1),
        (4.19, 4), (3.63, 4), (3.21, 2), (2.83, 2), (1.98, 2)
    ]
    
    print("Step 1: Analyzing the provided 1H NMR data.")
    
    # Calculate total integration and categorize protons
    total_protons_nmr = sum(d[1] for d in nmr_data)
    aromatic_protons_nmr = sum(d[1] for d in nmr_data if d[0] > 6.5)
    piperazine_protons_nmr = sum(d[1] for d in nmr_data if 3.5 <= d[0] <= 4.5)
    thq_aliphatic_protons_nmr = sum(d[1] for d in nmr_data if 1.5 <= d[0] < 3.5)

    print(f"Total observed protons from NMR integration: {total_protons_nmr}H")
    print("Breakdown of observed protons by type:")
    print(f"- Aromatic Protons (>6.5 ppm): {aromatic_protons_nmr}H")
    print(f"- Piperazine Protons (3.5-4.5 ppm): {piperazine_protons_nmr}H")
    print(f"- Tetrahydroquinoline Aliphatic Protons (1.5-3.5 ppm): {thq_aliphatic_protons_nmr}H")
    print("-" * 40)

    # 2. Define the expected proton counts for each compound
    compounds = {
        'A': {'total_H': 22, 'aromatic_H': 7, 'piperazine_H': 8, 'thq_aliphatic_H': 6, 'nh_H': 1, 'type': 'Ligand'},
        'C': {'total_H': 23, 'aromatic_H': 8, 'piperazine_H': 8, 'thq_aliphatic_H': 6, 'nh_H': 1, 'type': 'Ligand'},
        'B': {'total_H': 44, 'type': 'Complex [Zn(A)2]'},
        'D': {'total_H': 46, 'type': 'Complex [Zn(C)2]'},
        'E': {'total_H': 44, 'type': 'Complex [Zn(A)2]'} # Visually identical to B
    }
    
    print("Step 2: Calculating expected proton counts for each compound.")
    for name, data in compounds.items():
        print(f"Compound {name} ({data['type']}): Total Protons = {data['total_H']}", end="")
        if data['type'] == 'Ligand':
            print(f" (Aromatic={data['aromatic_H']}, NH={data['nh_H']})")
        else:
            print("")
    print("-" * 40)
    
    # 3. Compare NMR data with each compound and find the match
    print("Step 3: Comparing NMR data with each candidate compound.")
    match = None
    for name, data in compounds.items():
        # The total observed proton count (21) is much lower than the complexes (>40)
        if data['type'] != 'Ligand':
            continue
        
        # For ligands, check if observed protons match theoretical minus the exchangeable NH proton
        if total_protons_nmr == data['total_H'] - data['nh_H']:
            # Check the aromatic proton count, which is a key differentiator
            if aromatic_protons_nmr == data['aromatic_H']:
                # Check the other regions for confirmation
                if (piperazine_protons_nmr == data['piperazine_H'] and
                    thq_aliphatic_protons_nmr == data['thq_aliphatic_H']):
                    match = name
                    print(f"-> Match Found: Compound {name}")
                    print(f"   - The observed total of {total_protons_nmr}H matches Compound {name}'s count of {data['total_H']}H, assuming the single N-H proton is not observed.")
                    print(f"   - The observed aromatic count of {aromatic_protons_nmr}H perfectly matches the {data['aromatic_H']} aromatic protons of Compound {name}.")
                    print(f"   - The remaining proton distributions also match.")
                    break
    
    print("-" * 40)

    # 4. Final Conclusion
    print("Step 4: Final Conclusion.")
    if match:
        answer_choices = {'A': 'B', 'B': 'D', 'C': 'E', 'D': 'C', 'E': 'A'}
        final_answer_choice = [k for k, v in answer_choices.items() if v == match][0]
        
        print(f"The NMR data corresponds to the structure of Compound {match}.")
        print(f"In the multiple-choice list, Compound {match} is given as answer choice '{final_answer_choice}'.")
        
        # Display the data in an equation-like format
        print("\nSummary Equation of NMR Evidence:")
        equation_terms = [f"{d[0]}({d[1]}H)" for d in nmr_data]
        print(" + ".join(equation_terms) + f" = {total_protons_nmr} total protons")

    else:
        print("No definitive match could be found based on the analysis.")

identify_compound_from_nmr()