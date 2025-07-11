def solve_molecular_formula():
    """
    Calculates the molecular formula of product A based on the reaction described.
    The overall reaction replaces the -CH(COOMe)- fragment of compound 1 with a -CH(benzyl)- fragment.
    """

    # Molecular formula of Compound 1 (methyl 1-oxoindan-2-carboxylate) is C11H10O3
    formula_c1 = {'C': 11, 'H': 10, 'O': 3}

    # Formula of the fragment to be removed: -CH(COOMe)-
    # This consists of the C2 atom, its attached H, and the COOMe group.
    # CH (1C, 1H) + COOMe (2C, 3H, 2O)
    fragment_removed = {'C': 3, 'H': 4, 'O': 2}

    # Formula of the fragment to be added: -CH(benzyl)-
    # This consists of the C2 atom, its attached H, and the benzyl group.
    # CH (1C, 1H) + Benzyl (C7H7)
    fragment_added = {'C': 8, 'H': 8, 'O': 0}

    # Calculate the formula of the final product A
    final_C = formula_c1['C'] - fragment_removed['C'] + fragment_added['C']
    final_H = formula_c1['H'] - fragment_removed['H'] + fragment_added['H']
    final_O = formula_c1['O'] - fragment_removed['O'] + fragment_added['O']

    final_formula = f"C{final_C}H{final_H}O{final_O if final_O > 1 else ''}"

    print("To find the molecular formula of compound A, we calculate the net change from compound 1.")
    print("The net reaction replaces the -CH(COOMe)- fragment with a -CH(benzyl)- fragment.\n")

    print(f"Formula of Compound 1: C{formula_c1['C']}H{formula_c1['H']}O{formula_c1['O']}")
    print(f"Formula of removed fragment [-CH(COOMe)-]: C{fragment_removed['C']}H{fragment_removed['H']}O{fragment_removed['O']}")
    print(f"Formula of added fragment [-CH(benzyl)-]: C{fragment_added['C']}H{fragment_added['H']}\n")

    print("Calculating the number of atoms for each element in Product A:")
    print(f"Number of Carbon atoms   = {formula_c1['C']} - {fragment_removed['C']} + {fragment_added['C']} = {final_C}")
    print(f"Number of Hydrogen atoms = {formula_c1['H']} - {fragment_removed['H']} + {fragment_added['H']} = {final_H}")
    print(f"Number of Oxygen atoms   = {formula_c1['O']} - {fragment_removed['O']} + {fragment_added['O']} = {final_O}\n")

    print(f"The molecular formula of compound A is {final_formula}.")

solve_molecular_formula()