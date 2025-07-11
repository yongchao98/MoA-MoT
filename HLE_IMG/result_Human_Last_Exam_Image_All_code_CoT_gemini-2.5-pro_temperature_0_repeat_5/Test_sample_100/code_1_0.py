def solve_molecular_formula():
    """
    This script determines the molecular formula of the product of a three-step reaction.
    """

    # Step 0: Define the molecular formula of the starting material.
    # The starting material is 5-(trifluoromethyl)-2-(p-methoxybenzyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
    # Its formula is C15H14F3NO2.
    atoms_start = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    
    print("Step-by-step determination of the product's molecular formula:")
    print(f"Starting material: C{atoms_start['C']}H{atoms_start['H']}F{atoms_start['F']}NO{atoms_start['O']}")
    print("-" * 20)

    # Step 1: Reaction with CAN. This removes the p-methoxybenzyl (PMB) group (C8H9O) 
    # and replaces it with a hydrogen atom.
    # Change in atoms: -C8, -H9, -O1, +H1  =>  -C8, -H8, -O1
    c1 = atoms_start['C'] - 8
    h1 = atoms_start['H'] - 8
    o1 = atoms_start['O'] - 1
    print("Step 1: PMB deprotection removes C8H9O and adds H.")
    print(f"Equation: C{atoms_start['C']}H{atoms_start['H']}F{atoms_start['F']}NO{atoms_start['O']} -> C{c1}H{h1}F{atoms_start['F']}NO{o1}")
    print("-" * 20)

    # Step 2: Reaction with Pd/C, H2. This hydrogenation adds two hydrogen atoms.
    # Change in atoms: +H2
    h2 = h1 + 2
    print("Step 2: Hydrogenation adds H2.")
    print(f"Equation: C{c1}H{h1}F{atoms_start['F']}NO{o1} + H2 -> C{c1}H{h2}F{atoms_start['F']}NO{o1}")
    print("-" * 20)

    # Step 3: Reaction with 4 N HCl. This is a lactam hydrolysis, which consumes one molecule of water.
    # Change in atoms: +H2, +O1
    h_final = h2 + 2
    o_final = o1 + 1
    c_final = c1
    f_final = atoms_start['F']
    n_final = atoms_start['N']
    print("Step 3: Lactam hydrolysis adds H2O.")
    print(f"Equation: C{c1}H{h2}F{atoms_start['F']}NO{o1} + H2O -> C{c_final}H{h_final}F{f_final}NO{o_final}")
    print("-" * 20)

    # Construct the final formula string in standard chemical notation.
    final_formula_str = f"C{c_final}H{h_final}F{f_final}NO{o_final}"
    
    print(f"The final molecular formula of the product is: {final_formula_str}")

solve_molecular_formula()