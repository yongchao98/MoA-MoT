def solve_molecular_formula():
    """
    Calculates the molecular formula of product A by tracking atom changes
    through the multi-step reaction.
    """
    # 1. Define the molecular formula of Reactant 1:
    # Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate
    reactant1 = {'C': 11, 'H': 10, 'O': 3}

    # The group added from Reactant 2 (benzyl bromide) is a benzyl group.
    benzyl_group = {'C': 7, 'H': 7}

    # 2. Step 1: Alkylation
    # The acidic alpha-hydrogen is replaced by the benzyl group.
    # Net change: -H, +C7H7
    alkylation_product = reactant1.copy()
    alkylation_product['H'] -= 1  # Remove acidic H
    alkylation_product['C'] += benzyl_group['C']  # Add benzyl C
    alkylation_product['H'] += benzyl_group['H']  # Add benzyl H
    # Formula of the alkylated ester: C18H16O3

    # 3. Step 2: Saponification
    # The methyl ester (-COOCH3) is hydrolyzed to a carboxylic acid (-COOH).
    # Net change: replace -CH3 with -H.
    saponification_product = alkylation_product.copy()
    saponification_product['C'] -= 1  # Remove methyl C
    saponification_product['H'] -= 3  # Remove methyl H
    saponification_product['H'] += 1  # Add acid H
    # Formula of the beta-keto acid: C17H14O3

    # 4. Step 3: Decarboxylation
    # The beta-keto acid loses CO2.
    final_product_A = saponification_product.copy()
    final_product_A['C'] -= 1  # Remove C from CO2
    final_product_A['O'] -= 2  # Remove O from CO2
    # Final formula: C16H14O

    # 5. Print the result
    c = final_product_A['C']
    h = final_product_A['H']
    o = final_product_A['O']

    print("The molecular formula of compound A is calculated as follows:")
    print(f"Start with Reactant 1: C{reactant1['C']}H{reactant1['H']}O{reactant1['O']}")
    print(f"After alkylation: C{alkylation_product['C']}H{alkylation_product['H']}O{alkylation_product['O']}")
    print(f"After saponification: C{saponification_product['C']}H{saponification_product['H']}O{saponification_product['O']}")
    print(f"After decarboxylation: C{final_product_A['C']}H{final_product_A['H']}O{final_product_A['O']}")
    print("\nThe final molecular formula of compound A is:")
    print(f"Carbon atoms: {c}")
    print(f"Hydrogen atoms: {h}")
    print(f"Oxygen atoms: {o}")

solve_molecular_formula()
<<<C16H14O>>>