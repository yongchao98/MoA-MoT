def calculate_molecular_formula():
    """
    Calculates the molecular formula of product A based on the given reaction.
    """
    # Step 1: Define the molecular formula of the reactants.
    # Compound 1: Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate
    compound1 = {'C': 11, 'H': 10, 'O': 3}
    
    # The reaction adds a benzyl group (C7H7) from benzyl bromide.
    benzyl_group = {'C': 7, 'H': 7}

    print("Step 1: Alkylation")
    print("The reaction starts with compound 1 (C11H10O3).")
    print("A benzyl group (C7H7) is added and one hydrogen atom is removed.")
    
    # Calculate the formula of the alkylated intermediate
    intermediate_C = compound1['C'] + benzyl_group['C']
    intermediate_H = compound1['H'] - 1 + benzyl_group['H']
    intermediate_O = compound1['O']
    
    print(f"Carbon atoms in intermediate: {compound1['C']} + {benzyl_group['C']} = {intermediate_C}")
    print(f"Hydrogen atoms in intermediate: {compound1['H']} - 1 + {benzyl_group['H']} = {intermediate_H}")
    print(f"The formula of the alkylated ester intermediate is C{intermediate_C}H{intermediate_H}O{intermediate_O}.")
    print("-" * 20)

    print("Step 2: Saponification")
    print("The methyl ester (-COOCH3) is hydrolyzed to a carboxylic acid (-COOH).")
    print("This corresponds to replacing a methyl group (-CH3) with a hydrogen atom (-H), a net removal of CH2.")
    
    # Calculate the final product formula by removing CH2
    final_C = intermediate_C - 1
    final_H = intermediate_H - 2
    final_O = intermediate_O
    
    print(f"Carbon atoms in final product A: {intermediate_C} - 1 = {final_C}")
    print(f"Hydrogen atoms in final product A: {intermediate_H} - 2 = {final_H}")
    print(f"Oxygen atoms remain unchanged: {final_O}")
    print("-" * 20)
    
    # Display the final molecular formula
    final_formula = f"C{final_C}H{final_H}O{final_O}"
    print(f"The molecular formula of compound A is {final_formula}.")

calculate_molecular_formula()