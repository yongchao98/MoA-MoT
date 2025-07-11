def solve_molecular_formula():
    """
    This script calculates the molecular formula of Compound B based on the provided reaction scheme.
    """
    # Step 1: Define the atomic composition of the reactant cation (Cation-X).
    # From the structure (tetramethoxy-xanthylium cation), its formula is C17H17O5.
    cation_X = {'C': 17, 'H': 17, 'N': 0, 'O': 5}

    # Step 2: Define the atomic composition of the amine reagent (methyl-3-aminopropionate).
    # Structure H2N-CH2-CH2-COOCH3 gives the formula C4H9NO2.
    amine = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

    # Step 3: Define the atomic composition of the eliminated water molecule.
    water = {'C': 0, 'H': 2, 'N': 0, 'O': 1}

    # Step 4: Calculate the composition of Compound B's cation by atom balance.
    # The reaction is: Cation-X + Amine -> Cation-B + Water
    # So, Cation-B's atoms = (Cation-X's atoms + Amine's atoms) - Water's atoms.
    cation_B = {
        'C': cation_X['C'] + amine['C'] - water['C'],
        'H': cation_X['H'] + amine['H'] - water['H'],
        'N': cation_X['N'] + amine['N'] - water['N'],
        'O': cation_X['O'] + amine['O'] - water['O']
    }

    # Step 5: Format the results and print the full chemical equation.
    # The convention is to omit the number 1 for atoms with a count of one.
    cation_X_formula = f"C{cation_X['C']}H{cation_X['H']}O{cation_X['O']}"
    amine_formula = f"C{amine['C']}H{amine['H']}NO{amine['O']}"
    cation_B_formula = f"C{cation_B['C']}H{cation_B['H']}NO{cation_B['O']}"
    water_formula = f"H{water['H']}O"
    
    print("The overall reaction to form the cation of Compound B is:")
    print(f"{cation_X_formula} + {amine_formula} -> {cation_B_formula} + {water_formula}")

    print("\nTherefore, the molecular formula of the cation of compound B is:")
    print(cation_B_formula)

solve_molecular_formula()