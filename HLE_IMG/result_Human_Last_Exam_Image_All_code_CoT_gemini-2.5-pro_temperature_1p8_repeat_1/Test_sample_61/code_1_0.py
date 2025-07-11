def calculate_molecular_formula():
    """
    Calculates the molecular formula of compound A based on the reaction provided.
    """

    # Molecular formula of Compound 1 (methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate)
    # C=11, H=10, O=3
    compound1 = {'C': 11, 'H': 10, 'O': 3}

    # Atomic composition of the groups involved in the reaction
    # This is the benzyl group from compound 2 (C7H7Br) that is added.
    benzyl_group = {'C': 7, 'H': 7}
    # This represents the methoxycarbonyl group (-COOCH3) that is ultimately lost.
    methoxycarbonyl_group = {'C': 2, 'H': 3, 'O': 2}
    # An alpha-hydrogen is replaced, and another hydrogen is gained in its place
    # after decarboxylation, so the net change to H from this is 0.

    # Calculate the final formula of product A.
    # The reaction can be summarized as:
    # Compound1 + Benzyl_group - Methoxycarbonyl_group = Compound A
    final_formula = {}
    final_formula['C'] = compound1['C'] + benzyl_group['C'] - methoxycarbonyl_group['C']
    # The hydrogen calculation: H starts at 10. Lose 1 alpha-H, gain 7 benzyl-H, lose 3 methyl-H from ester, gain 1 H from hydrolysis/workup.
    # H_final = 10 - 1 + 7 - 3 + 1 = 14
    final_formula['H'] = compound1['H'] + benzyl_group['H'] - methoxycarbonyl_group['H']
    final_formula['O'] = compound1['O'] - methoxycarbonyl_group['O']

    c = final_formula['C']
    h = final_formula['H']
    o = final_formula['O']

    # The prompt asks to output each number in the final equation.
    # Here, we construct the final molecular formula string and print its components.
    print("The final molecular formula of compound A is:")
    print("C", c, "H", h, "O", o, sep="")


calculate_molecular_formula()