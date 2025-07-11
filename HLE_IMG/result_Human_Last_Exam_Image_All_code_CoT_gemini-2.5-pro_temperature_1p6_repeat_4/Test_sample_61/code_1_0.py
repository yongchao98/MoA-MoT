def solve_molecular_formula():
    """
    Calculates the molecular formula of product A based on the reaction provided.
    The reaction is an alkylation followed by saponification and decarboxylation.
    """
    # Step 1: Define the molecular formulas of reactants and relevant groups as dictionaries.
    # Compound 1 is methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate: C11H10O3
    compound_1 = {'C': 11, 'H': 10, 'O': 3}
    
    # The reaction adds a benzyl group (C7H7) from benzyl bromide
    # and removes a hydrogen atom in the alkylation step.
    benzyl_group = {'C': 7, 'H': 7}
    hydrogen_atom = {'H': 1}
    
    # The reaction also involves replacing the carbo-methoxy group (-COOMe) with a hydrogen atom.
    # COOMe group: C2H3O2
    coome_group = {'C': 2, 'H': 3, 'O': 2}

    print("Step-by-step calculation of the molecular formula for product A:")
    print("-" * 60)
    print(f"1. The molecular formula of reactant 1 is C{compound_1['C']}H{compound_1['H']}O{compound_1['O']}.")

    # Step 2: Calculate the formula of the alkylated intermediate.
    # This is (Compound 1) - (Hydrogen) + (Benzyl group)
    intermediate_C = compound_1['C'] + benzyl_group['C']
    intermediate_H = compound_1['H'] - hydrogen_atom['H'] + benzyl_group['H']
    intermediate_O = compound_1['O']
    
    print("\n2. First, an alkylation reaction occurs. A hydrogen atom is replaced by a benzyl group (C7H7).")
    print(f"   Formula of intermediate = C({compound_1['C']})H({compound_1['H']})O({compound_1['O']}) - H + C({benzyl_group['C']})H({benzyl_group['H']})")
    print(f"                         = C({intermediate_C})H({intermediate_H})O({intermediate_O})")

    # Step 3: Calculate the formula of the final product A.
    # This is (Intermediate) - (COOMe group) + (Hydrogen)
    product_A_C = intermediate_C - coome_group['C']
    product_A_H = intermediate_H - coome_group['H'] + hydrogen_atom['H']
    product_A_O = intermediate_O - coome_group['O']

    print("\n3. Next, saponification and decarboxylation occur, which replaces the -COOMe group (C2H3O2) with a hydrogen atom.")
    print(f"   Formula of Product A = C({intermediate_C})H({intermediate_H})O({intermediate_O}) - C({coome_group['C']})H({coome_group['H']})O({coome_group['O']}) + H")
    # To avoid printing O1, we use a small helper function
    def format_formula(atoms):
        c, h, o = atoms['C'], atoms['H'], atoms['O']
        o_str = f"O({o})" if o != 0 else ""
        return f"C({c})H({h}){o_str}"

    print(f"                        = C({product_A_C})H({product_A_H})O({product_A_O})")
    
    # Step 4: Format the final formula string
    o_str_final = f"O" if product_A_O == 1 else f"O{product_A_O}" if product_A_O > 1 else ""
    final_formula = f"C{product_A_C}H{product_A_H}{o_str_final}"
    
    print("-" * 60)
    print(f"The final molecular formula of compound A is {final_formula}.")

solve_molecular_formula()