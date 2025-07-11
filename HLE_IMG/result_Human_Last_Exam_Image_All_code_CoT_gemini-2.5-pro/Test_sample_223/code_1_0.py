def calculate_molecular_formula():
    """
    This script calculates the molecular formula of compound B based on the provided reaction scheme.
    """
    
    # Step 1: Define the atomic composition of the reactants and byproducts.
    # The starting material is 9-(2,6-dimethoxyphenyl)-1,8-dimethoxyxanthylium cation.
    start_cation = {'C': 23, 'H': 21, 'O': 5, 'N': 0}
    
    # The reagent is methyl-3-aminopropionate (NH2-CH2-CH2-COOCH3).
    amine = {'C': 4, 'H': 9, 'N': 1, 'O': 2}
    
    # The reaction is a condensation that eliminates one molecule of water.
    eliminated_molecule = {'C': 0, 'H': 2, 'O': 1, 'N': 0}
    
    # Print the initial information
    print("Step-by-step calculation of the molecular formula for the cation of Compound B:")
    print("-" * 70)
    print(f"1. Formula of Starting Cation: C{start_cation['C']}H{start_cation['H']}O{start_cation['O']}")
    print(f"2. Formula of Amine Reagent (methyl-3-aminopropionate): C{amine['C']}H{amine['H']}N{amine['N']}O{amine['O']}")
    print(f"3. Formula of Eliminated Molecule (water): H{eliminated_molecule['H']}O{eliminated_molecule['O']}")
    print("-" * 70)
    
    # Step 2: Calculate the atomic composition of the product cation.
    # Formula = (Start Cation) + (Amine) - (Water)
    product_cation = {}
    product_cation['C'] = start_cation['C'] + amine['C'] - eliminated_molecule['C']
    product_cation['H'] = start_cation['H'] + amine['H'] - eliminated_molecule['H']
    product_cation['N'] = start_cation['N'] + amine['N'] - eliminated_molecule['N']
    product_cation['O'] = start_cation['O'] + amine['O'] - eliminated_molecule['O']
    
    # Print the calculation steps
    print("4. Calculation for each element in the product cation:")
    print(f"   Carbon (C)   = {start_cation['C']} + {amine['C']} - {eliminated_molecule['C']} = {product_cation['C']}")
    print(f"   Hydrogen (H) = {start_cation['H']} + {amine['H']} - {eliminated_molecule['H']} = {product_cation['H']}")
    print(f"   Nitrogen (N) = {start_cation['N']} + {amine['N']} - {eliminated_molecule['N']} = {product_cation['N']}")
    print(f"   Oxygen (O)   = {start_cation['O']} + {amine['O']} - {eliminated_molecule['O']} = {product_cation['O']}")
    print("-" * 70)
    
    # Step 3: Construct and print the final molecular formula.
    # The organic product is an acridinium cation.
    cation_formula_str = f"C{product_cation['C']}H{product_cation['H']}N{product_cation['N']}O{product_cation['O']}"
    print(f"The molecular formula of the organic cation in Compound B is: {cation_formula_str}")
    
    # Note on the full compound: Compound B is the salt with the BF4- counterion.
    # The full formula is C27H28BF4NO6. The question is interpreted as asking for the organic part.
    print("\nNote: Compound B is the salt with the tetrafluoroborate (BF4-) anion.")
    print(f"The full molecular formula of Compound B is C{product_cation['C']}H{product_cation['H']}BF4N{product_cation['N']}O{product_cation['O']}.")

if __name__ == '__main__':
    calculate_molecular_formula()