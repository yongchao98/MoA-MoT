def solve_molecular_formula():
    """
    This function calculates the molecular formula of the final product based on the provided reaction scheme.
    """
    # Step 1: Determine the molecular formula of the starting material.
    # Structure: 2-(p-methoxybenzyl)-6-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
    # Bicyclic core (C6H5F3NO): C=6, H=5, F=3, N=1, O=1 (from C=O).
    # PMB group (-C8H9O): C=8, H=9, O=1.
    # Total atoms in Starting Material: C=14, H=14, F=3, N=1, O=2.
    starting_material = {'C': 14, 'H': 14, 'F': 3, 'N': 1, 'O': 2}

    # Step 2: Calculate the formula for Intermediate 1.
    # Reaction: Deprotection of PMB group using CAN. This removes the C8H9O group and adds an H atom.
    # Let's calculate the net change in atoms: remove PMB group (-C8H9), and its oxygen (-O), and add one hydrogen (-H).
    # This leads to a net change of -C8, H(-9+1)=-8, -O1. So, -C8H8O.
    # A simpler way is to count atoms of the deprotected structure: 6-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
    # Atoms: C=6, H=6, F=3, N=1, O=1.
    intermediate_1 = {'C': 6, 'H': 6, 'F': 3, 'N': 1, 'O': 1}

    # Step 3: Calculate the formula for Intermediate 2.
    # Reaction: Hydrogenation with Pd/C and H2. This reduces the C=C double bond.
    # It adds two hydrogen atoms (H2) to Intermediate 1.
    intermediate_2 = intermediate_1.copy()
    intermediate_2['H'] += 2
    # Formula of Intermediate 2: C6H8F3NO.

    # Step 4: Calculate the formula for the final Product.
    # Reaction: Hydrolysis of the lactam with hot acid.
    # This reaction opens the cyclic amide ring, which corresponds to the net addition of one water molecule (H2O).
    product = intermediate_2.copy()
    product['H'] += 2
    product['O'] += 1
    # Final Formula: C6H10F3NO2.

    # Step 5: Print the final molecular formula and the count of each atom.
    c = product['C']
    h = product['H']
    f = product['F']
    n = product['N']
    o = product['O']

    formula_string = f"C{c}H{h}F{f}NO{o}" # N is 1, so no number needed

    print(f"The molecular formula of the final product is {formula_string}")
    print("The individual atom counts in the final product are:")
    print(f"Carbon (C): {c}")
    print(f"Hydrogen (H): {h}")
    print(f"Fluorine (F): {f}")
    print(f"Nitrogen (N): {n}")
    print(f"Oxygen (O): {o}")

solve_molecular_formula()