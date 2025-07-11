def get_molecular_formula():
    """
    This script calculates the molecular formula of the final product A
    by breaking down the reaction into three steps: alkylation, saponification,
    and decarboxylation, and tracking the atomic composition at each stage.
    """
    # Step 1: Define the molecular formula of the starting materials.
    # Compound 1 is Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate.
    # Its formula is C11H10O3.
    c1, h1, o1 = 11, 10, 3

    # Compound 2 is Benzyl bromide.
    # Its formula is C7H7Br.
    c2, h2 = 7, 7

    print("--- Reaction Analysis ---")
    print(f"Molecular formula of Compound 1: C{c1}H{h1}O{o1}")
    print(f"Molecular formula of Compound 2: C{c2}H{h2}Br")
    print("-" * 30)

    # Step 2: Calculate the formula of the alkylation product.
    # The reaction replaces the acidic alpha-hydrogen of Compound 1 with the benzyl group (C7H7) from Compound 2.
    c_alk = c1 + c2
    h_alk = h1 - 1 + h2
    o_alk = o1
    print("Step A: Alkylation")
    print(f"Formula after alkylation = (C{c1}H{h1}O{o1} - H) + (C{c2}H{h2})")
    print(f"Equation: C({c1}+{c2}) H({h1}-1+{h2}) O({o1}) = C{c_alk}H{h_alk}O{o_alk}")
    print("-" * 30)

    # Step 3: Calculate the formula after saponification.
    # The methyl ester group (-COOMe) is hydrolyzed to a carboxylic acid (-COOH).
    # This corresponds to a net change of replacing -CH3 with -H, or a loss of CH2.
    c_acid = c_alk - 1
    h_acid = h_alk - 2
    o_acid = o_alk
    print("Step B: Saponification (Hydrolysis)")
    print(f"The formula of the resulting β-keto acid is C{c_acid}H{h_acid}O{o_acid}.")
    print(f"Equation: C({c_alk}-1) H({h_alk}-2) O({o_alk}) = C{c_acid}H{h_acid}O{o_acid}")
    print("-" * 30)

    # Step 4: Calculate the formula of the final product A after decarboxylation.
    # The β-keto acid loses a molecule of carbon dioxide (CO2).
    c_A = c_acid - 1
    h_A = h_acid
    o_A = o_acid - 2
    print("Step C: Decarboxylation")
    print(f"The final product A is formed by losing CO2 from the β-keto acid.")
    print(f"Equation: C({c_acid}-1) H({h_acid}) O({o_acid}-2) = C{c_A}H{h_A}O{o_A}")
    print("-" * 30)

    # Step 5: Display the final result.
    # Check if oxygen is 1, to omit the subscript.
    if o_A == 1:
        final_formula = f"C{c_A}H{h_A}O"
    else:
        final_formula = f"C{c_A}H{h_A}O{o_A}"

    print(f"The final molecular formula of compound A is {final_formula}.")


# Execute the function to find the answer.
get_molecular_formula()