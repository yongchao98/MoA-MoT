def calculate_product_formula():
    """
    Calculates the molecular formula of product A based on the reaction stoichiometry.
    The reaction is an alkylation of a beta-keto ester followed by saponification
    and decarboxylation.
    """

    # Molecular formula of Compound 1 (C11H10O3)
    compound1_atoms = {'C': 11, 'H': 10, 'O': 3}

    # Atoms from the benzyl group (C7H7) that is added
    benzyl_group_atoms = {'C': 7, 'H': 7}

    # Atoms from the methoxycarbonyl group (COOMe) that is removed (C2H3O2)
    coome_group_atoms = {'C': 2, 'H': 3, 'O': 2}

    # The overall transformation is the replacement of the alpha-H and the -COOMe group
    # at C2 with a benzyl group and a new H atom.
    # Net change: -H(alpha) + C7H7(benzyl) - C2H3O2(COOMe) + H(final)

    # Calculate the final number of atoms for compound A
    # Initial atoms from Compound 1
    c_final = compound1_atoms['C']
    h_final = compound1_atoms['H']
    o_final = compound1_atoms['O']
    
    # Equation for Carbon atoms
    c_eq_str = f"x (Carbon) = {c_final} + {benzyl_group_atoms['C']} (from benzyl) - {coome_group_atoms['C']} (from COOMe) = {c_final + benzyl_group_atoms['C'] - coome_group_atoms['C']}"
    c_final = c_final + benzyl_group_atoms['C'] - coome_group_atoms['C']

    # Equation for Hydrogen atoms
    h_eq_str = f"y (Hydrogen) = {h_final} - 1 (alpha-H) + {benzyl_group_atoms['H']} (from benzyl) - {coome_group_atoms['H']} (from COOMe) + 1 (final-H) = {h_final - 1 + benzyl_group_atoms['H'] - coome_group_atoms['H'] + 1}"
    h_final = h_final - 1 + benzyl_group_atoms['H'] - coome_group_atoms['H'] + 1
    
    # Equation for Oxygen atoms
    o_eq_str = f"z (Oxygen) = {o_final} - {coome_group_atoms['O']} (from COOMe) = {o_final - coome_group_atoms['O']}"
    o_final = o_final - coome_group_atoms['O']

    # Print the step-by-step calculation for each element
    print("Calculation for the molecular formula of compound A (CxHyOz):")
    print(c_eq_str)
    print(h_eq_str)
    print(o_eq_str)

    # Print the final molecular formula
    print(f"\nThe molecular formula of compound A is C{c_final}H{h_final}O{o_final}.")

calculate_product_formula()