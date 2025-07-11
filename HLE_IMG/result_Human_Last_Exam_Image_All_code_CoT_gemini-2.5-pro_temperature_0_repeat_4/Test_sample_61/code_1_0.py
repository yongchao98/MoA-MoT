def calculate_product_formula():
    """
    Calculates the molecular formula of product A based on the reaction sequence.
    """
    # Step 0: Define initial molecular formulas and groups.
    # Compound 1: Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate
    # C (benzene ring) = 6, C (5-membered ring) = 3, C (ester) = 2 -> Total C = 11
    # H (benzene ring) = 4, H (C2) = 1, H (C3) = 2, H (methyl) = 3 -> Total H = 10
    # O (ketone) = 1, O (ester) = 2 -> Total O = 3
    # Formula: C11H10O3
    compound1 = {'C': 11, 'H': 10, 'O': 3}

    # Benzyl group (C6H5CH2-)
    benzyl_group = {'C': 7, 'H': 7}
    
    # Hydrogen atom
    hydrogen_atom = {'H': 1}

    # Methyl group (-CH3)
    methyl_group = {'C': 1, 'H': 3}

    # Carbon dioxide molecule (CO2)
    co2_molecule = {'C': 1, 'O': 2}

    print("Starting with Compound 1:")
    print(f"C{compound1['C']} + H{compound1['H']} + O{compound1['O']}")
    print("-" * 30)

    # Step 1: Alkylation
    # Replace one acidic H atom with a benzyl group.
    # Intermediate = Compound1 - H + Benzyl
    alkylated_intermediate = compound1.copy()
    alkylated_intermediate['H'] -= hydrogen_atom['H']
    alkylated_intermediate['C'] += benzyl_group['C']
    alkylated_intermediate['H'] += benzyl_group['H']
    print("Step 1: Alkylation (replaces H with Benzyl group)")
    print(f"Formula of alkylated intermediate: C{alkylated_intermediate['C']}H{alkylated_intermediate['H']}O{alkylated_intermediate['O']}")
    print(f"Calculation: C({compound1['C']}+{benzyl_group['C']}) H({compound1['H']}-{hydrogen_atom['H']}+{benzyl_group['H']}) O({compound1['O']})")
    print("-" * 30)

    # Step 2: Saponification (Ester Hydrolysis)
    # The methyl ester (-COOMe) becomes a carboxylic acid (-COOH).
    # This is equivalent to replacing a methyl group (-CH3) with a hydrogen atom (-H).
    # Intermediate = Alkylated Intermediate - CH3 + H
    beta_keto_acid = alkylated_intermediate.copy()
    beta_keto_acid['C'] -= methyl_group['C']
    beta_keto_acid['H'] -= methyl_group['H']
    beta_keto_acid['H'] += hydrogen_atom['H']
    print("Step 2: Saponification (replaces -CH3 with -H)")
    print(f"Formula of beta-keto acid: C{beta_keto_acid['C']}H{beta_keto_acid['H']}O{beta_keto_acid['O']}")
    print(f"Calculation: C({alkylated_intermediate['C']}-{methyl_group['C']}) H({alkylated_intermediate['H']}-{methyl_group['H']}+{hydrogen_atom['H']}) O({alkylated_intermediate['O']})")
    print("-" * 30)

    # Step 3: Decarboxylation
    # The beta-keto acid loses a molecule of CO2.
    # Product A = Beta-keto acid - CO2
    product_A = beta_keto_acid.copy()
    product_A['C'] -= co2_molecule['C']
    product_A['O'] -= co2_molecule['O']
    print("Step 3: Decarboxylation (loss of CO2)")
    print(f"Final formula of Compound A: C{product_A['C']}H{product_A['H']}O{product_A['O']}")
    print(f"Calculation: C({beta_keto_acid['C']}-{co2_molecule['C']}) H({beta_keto_acid['H']}) O({beta_keto_acid['O']}-{co2_molecule['O']})")
    print("-" * 30)

calculate_product_formula()