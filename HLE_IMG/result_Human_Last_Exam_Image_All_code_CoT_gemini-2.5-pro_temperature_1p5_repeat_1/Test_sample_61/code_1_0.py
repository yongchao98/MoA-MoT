def solve_chemistry_problem():
    """
    This function calculates the molecular formula of product A based on the reaction provided.
    The reaction involves alkylation, saponification, and decarboxylation.
    """
    
    # Molecular formulas of the reactants
    reactant1 = {'C': 11, 'H': 10, 'O': 3}
    reactant2 = {'C': 7, 'H': 7, 'Br': 1}
    
    # Groups involved in the transformation
    hbr = {'H': 1, 'Br': 1}
    coome_group = {'C': 2, 'H': 3, 'O': 2}
    h_atom = {'H': 1}
    
    print("Step 1: Determine the molecular formulas of the reactants.")
    print(f"Formula of reactant 1 (methyl 1-oxoindane-2-carboxylate): C{reactant1['C']}H{reactant1['H']}O{reactant1['O']}")
    print(f"Formula of reactant 2 (benzyl bromide): C{reactant2['C']}H{reactant2['H']}Br")
    print("-" * 30)
    
    print("Step 2: Calculate the formula of the alkylation intermediate.")
    print("This step combines reactant 1 and 2 and eliminates HBr.")
    
    inter_c = reactant1['C'] + reactant2['C']
    inter_h = reactant1['H'] + reactant2['H'] - hbr['H']
    inter_o = reactant1['O']
    
    print(f"Intermediate Carbon atoms: {reactant1['C']} + {reactant2['C']} = {inter_c}")
    print(f"Intermediate Hydrogen atoms: {reactant1['H']} + {reactant2['H']} - {hbr['H']} = {inter_h}")
    print(f"Intermediate Oxygen atoms: {reactant1['O']}")
    print(f"Intermediate formula: C{inter_c}H{inter_h}O{inter_o}")
    print("-" * 30)

    print("Step 3: Calculate the formula of the final product A.")
    print("This step involves saponification and decarboxylation, which removes the -COOMe group and adds an H atom.")
    
    final_c = inter_c - coome_group['C']
    final_h = inter_h - coome_group['H'] + h_atom['H']
    final_o = inter_o - coome_group['O']

    print(f"Final Carbon atoms: {inter_c} - {coome_group['C']} = {final_c}")
    print(f"Final Hydrogen atoms: {inter_h} - {coome_group['H']} + {h_atom['H']} = {final_h}")
    print(f"Final Oxygen atoms: {inter_o} - {coome_group['O']} = {final_o}")
    print("-" * 30)
    
    print(f"The final molecular formula of compound A is: C{final_c}H{final_h}O{final_o}")

solve_chemistry_problem()