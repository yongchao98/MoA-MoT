def calculate_molecular_formula():
    """
    Calculates the molecular formula of product A by tracking atomic changes
    through the reaction steps.
    """
    # Molecular formula of compound 1 (C11H10O3)
    c1 = {'C': 11, 'H': 10, 'O': 3}

    # Step 1: Alkylation
    # A proton (-H) is removed and a benzyl group (C7H7) is added.
    # Net change is addition of C7H6.
    alkylated_intermediate_C = c1['C'] + 7
    alkylated_intermediate_H = c1['H'] - 1 + 7
    alkylated_intermediate_O = c1['O']

    # Step 2: Saponification
    # The ester (-COOCH3) is hydrolyzed to an acid (-COOH).
    # This corresponds to replacing a methyl group (-CH3) with a hydrogen (-H).
    # Net change is loss of CH2.
    keto_acid_C = alkylated_intermediate_C - 1
    keto_acid_H = alkylated_intermediate_H - 2
    keto_acid_O = alkylated_intermediate_O

    # Step 3: Decarboxylation
    # The beta-keto acid loses a molecule of carbon dioxide (-CO2).
    final_product_C = keto_acid_C - 1
    final_product_H = keto_acid_H
    final_product_O = keto_acid_O - 2

    # Print the step-by-step calculation
    print("Derivation of the molecular formula for compound A:")
    print(f"1. Starting with Compound 1: C{c1['C']}H{c1['H']}O{c1['O']}")
    print(f"2. After Alkylation (adding C7H7, removing H):")
    print(f"   C = {c1['C']} + 7 = {alkylated_intermediate_C}")
    print(f"   H = {c1['H']} - 1 + 7 = {alkylated_intermediate_H}")
    print(f"   O = {c1['O']}")
    print(f"   Intermediate formula: C{alkylated_intermediate_C}H{alkylated_intermediate_H}O{alkylated_intermediate_O}")
    print(f"3. After Saponification (replacing -CH3 with -H):")
    print(f"   C = {alkylated_intermediate_C} - 1 = {keto_acid_C}")
    print(f"   H = {alkylated_intermediate_H} - 2 = {keto_acid_H}")
    print(f"   O = {alkylated_intermediate_O}")
    print(f"   Intermediate formula: C{keto_acid_C}H{keto_acid_H}O{keto_acid_O}")
    print(f"4. After Decarboxylation (removing CO2):")
    print(f"   C = {keto_acid_C} - 1 = {final_product_C}")
    print(f"   H = {keto_acid_H} = {final_product_H}")
    print(f"   O = {keto_acid_O} - 2 = {final_product_O}")
    print("-" * 20)
    print(f"Final molecular formula of compound A: C{final_product_C}H{final_product_H}O{final_product_O}")


calculate_molecular_formula()