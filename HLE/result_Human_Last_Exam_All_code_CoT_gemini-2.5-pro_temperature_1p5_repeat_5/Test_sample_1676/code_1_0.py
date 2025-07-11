def solve_chemistry_problem():
    """
    This function follows the reaction sequence to determine the molecular formula of the final product.
    """
    # The starting material is Terpinolene (p-menth-1,4(8)-diene).
    # Its molecular formula is C10H16.
    # We use a dictionary to keep track of the count of each atom.
    composition = {'C': 10, 'H': 16, 'O': 0, 'S': 0}
    print(f"Starting material: Terpinolene (Formula: C{composition['C']}H{composition['H']})")
    print("-" * 30)

    # Step 1: Reaction with m-CPBA (Epoxidation)
    # An oxygen atom is added to the molecule.
    composition['O'] += 1
    print("Step 1: Epoxidation with m-CPBA.")
    print("Reaction adds 1 Oxygen atom.")
    print(f"Compound 1 formula: C{composition['C']}H{composition['H']}O{composition['O']}")
    print("-" * 30)

    # Step 2: Reaction with N,N-dimethyl thioformamide
    # This converts the epoxide to a thiirane, replacing the Oxygen with a Sulfur atom.
    composition['O'] -= 1
    composition['S'] += 1
    print("Step 2: Thiirane formation with N,N-dimethyl thioformamide.")
    print("Reaction replaces 1 Oxygen with 1 Sulfur atom.")
    print(f"Compound 2 formula: C{composition['C']}H{composition['H']}S{composition['S']}")
    print("-" * 30)

    # Step 3: Reduction with LiAlH4
    # The thiirane is reductively opened, which adds two Hydrogen atoms to the molecule.
    composition['H'] += 2
    print("Step 3: Reduction of thiirane with LiAlH4.")
    print("Reaction adds 2 Hydrogen atoms.")
    print(f"Compound 3 formula: C{composition['C']}H{composition['H']}S{composition['S']}")
    print("-" * 30)
    
    # Final Result
    final_product_name = "p-menth-4(8)-en-1-thiol"
    formula_c = composition['C']
    formula_h = composition['H']
    formula_s = composition['S']
    
    print(f"The final product, Compound 3, is {final_product_name}.")
    print(f"The final molecular formula is C{formula_c}H{formula_h}S{formula_s}.")
    print("\nThe numbers in the final formula are:")
    print(f"Number of Carbon atoms: {formula_c}")
    print(f"Number of Hydrogen atoms: {formula_h}")
    print(f"Number of Sulfur atoms: {formula_s}")

solve_chemistry_problem()