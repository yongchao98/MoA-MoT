def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product by tracking atomic changes through the reaction scheme.
    """
    # Step 1: Determine the molecular formula of the starting material.
    # Bicyclic core (2-azabicyclo[2.2.1]hept-5-en-3-one): C6H5NO
    # Trifluoromethyl group (CF3): CF3
    # p-Methoxybenzyl (PMB) group (C6H4-CH2-OCH3): C8H9O
    atoms = {
        'C': 6 + 1 + 8,  # Core + CF3 + PMB
        'H': 5 + 9,      # Core + PMB
        'F': 3,          # CF3
        'N': 1,          # Core
        'O': 1 + 1       # Core + PMB
    }
    # print(f"Starting Material (C{atoms['C']}H{atoms['H']}F{atoms['F']}N{atoms['N']}O{atoms['O']}):")

    # Step 2: Reaction 1 -> Intermediate 1 (Deprotection of PMB)
    # Remove PMB group (C8H9O) and add 1 H to the nitrogen.
    atoms['C'] -= 8
    atoms['H'] -= 9
    atoms['O'] -= 1
    atoms['H'] += 1
    # print(f"Intermediate 1 (C{atoms['C']}H{atoms['H']}F{atoms['F']}N{atoms['N']}O{atoms['O']}):")

    # Step 3: Reaction 2 -> Intermediate 2 (Hydrogenation)
    # Add H2 to reduce the C=C double bond.
    atoms['H'] += 2
    # print(f"Intermediate 2 (C{atoms['C']}H{atoms['H']}F{atoms['F']}N{atoms['N']}O{atoms['O']}):")

    # Step 4: Reaction 3 -> Product (Lactam Hydrolysis)
    # Add one molecule of water (H2O).
    atoms['H'] += 2
    atoms['O'] += 1
    # print("Final Product:")

    # Final step: Print the number of atoms for each element in the final product.
    c_final = atoms['C']
    h_final = atoms['H']
    f_final = atoms['F']
    n_final = atoms['N']
    o_final = atoms['O']

    print(f"The final product contains:")
    print(f"Carbon (C) atoms: {c_final}")
    print(f"Hydrogen (H) atoms: {h_final}")
    print(f"Fluorine (F) atoms: {f_final}")
    print(f"Nitrogen (N) atoms: {n_final}")
    print(f"Oxygen (O) atoms: {o_final}")

    # Construct the molecular formula string
    formula = f"C{c_final}H{h_final}F{f_final}N{n_final}O{o_final}"
    # To handle elements with a count of 1 (e.g., N1 -> N)
    formula = formula.replace('F1', 'F').replace('N1', 'N').replace('O1', 'O')
    
    print(f"\nThe molecular formula of the product is: {formula}")

calculate_molecular_formula()
<<<C7H10F3NO2>>>