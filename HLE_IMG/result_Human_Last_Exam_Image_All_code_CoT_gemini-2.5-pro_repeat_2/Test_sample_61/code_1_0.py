def calculate_molecular_formula():
    """
    Calculates the molecular formula of product A based on the reaction sequence.
    """
    # Step 1: Define the molecular formula of the starting material, Compound 1.
    # Compound 1 is Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate.
    # C9H7O (1-oxoindan-2-yl skeleton) + C2H3O2 (methoxycarbonyl group) = C11H10O3
    compound1 = {'C': 11, 'H': 10, 'O': 3}
    print(f"The molecular formula of starting Compound 1 is C{compound1['C']}H{compound1['H']}O{compound1['O']}.")

    # Step 2: Define the groups being added or removed.
    # Alkylation adds a benzyl group (C7H7) and removes one proton (H).
    benzyl_group = {'C': 7, 'H': 7, 'O': 0}
    removed_proton = {'C': 0, 'H': 1, 'O': 0}
    
    # Saponification and decarboxylation removes the methoxycarbonyl group (-COOMe)
    # and adds one proton (H).
    coome_group = {'C': 2, 'H': 3, 'O': 2}
    added_proton = {'C': 0, 'H': 1, 'O': 0}

    print("\nThe reaction is an alkylation followed by saponification and decarboxylation.")
    print(f"This process effectively replaces the -COOMe group (C{coome_group['C']}H{coome_group['H']}O{coome_group['O']}) with a proton (H),")
    print(f"and the alpha-proton (H) with a benzyl group (C{benzyl_group['C']}H{benzyl_group['H']}).")
    
    # Step 3: Calculate the final molecular formula of product A.
    # We can combine the operations:
    # Final = Compound1 - removed_H + Benzyl - COOMe + added_H
    
    # Calculate the number of Carbon atoms in A
    final_C = compound1['C'] - removed_proton['C'] + benzyl_group['C'] - coome_group['C'] + added_proton['C']
    
    # Calculate the number of Hydrogen atoms in A
    final_H = compound1['H'] - removed_proton['H'] + benzyl_group['H'] - coome_group['H'] + added_proton['H']
    
    # Calculate the number of Oxygen atoms in A
    final_O = compound1['O'] - removed_proton['O'] + benzyl_group['O'] - coome_group['O'] + added_proton['O']

    print("\n--- Calculation of the Final Molecular Formula ---")
    print("Carbon atoms:")
    print(f"  {compound1['C']} (from Cpd 1) + {benzyl_group['C']} (from Benzyl) - {coome_group['C']} (from COOMe) = {final_C}")
    print("Hydrogen atoms:")
    print(f"  {compound1['H']} (from Cpd 1) - {removed_proton['H']} (removed H) + {benzyl_group['H']} (from Benzyl) - {coome_group['H']} (from COOMe) + {added_proton['H']} (added H) = {final_H}")
    print("Oxygen atoms:")
    print(f"  {compound1['O']} (from Cpd 1) - {coome_group['O']} (from COOMe) = {final_O}")

    print(f"\nThe final molecular formula of compound A is C{final_C}H{final_H}O{final_O}.")

calculate_molecular_formula()