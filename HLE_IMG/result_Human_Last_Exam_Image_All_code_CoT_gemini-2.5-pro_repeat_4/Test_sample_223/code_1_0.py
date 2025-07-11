def calculate_molecular_formula_B():
    """
    Calculates the molecular formula of compound B's organic cation based on the reaction scheme.
    """

    # Step 1: Define the atomic composition of the starting xanthylium cation (S+).
    # From structure analysis: 9-(2,6-dimethoxyphenyl)-1,8-dimethoxyxanthylium cation
    # Top group (C6H3(OCH3)2): C=8, H=9, O=2
    # Xanthylium part (C13H6O(OCH3)2+): C=15, H=12, O=3
    # Total S+ = C(8+15)H(9+12)O(2+3)+ = C23H21O5+
    s_cation = {'C': 23, 'H': 21, 'O': 5, 'N': 0}
    print(f"Step 1: The molecular formula of the starting cation S+ is C{s_cation['C']}H{s_cation['H']}O{s_cation['O']}.")

    # Step 2: Define the atomic composition of the amine reagent, methyl-3-aminopropionate.
    # Formula: H2N-CH2-CH2-COOCH3
    # C = 1(CH2) + 1(CH2) + 1(COO) + 1(CH3) = 4
    # H = 2(NH2) + 2(CH2) + 2(CH2) + 3(CH3) = 9
    # N = 1
    # O = 2(COO)
    amine = {'C': 4, 'H': 9, 'O': 2, 'N': 1}
    print(f"Step 2: The molecular formula of the reagent methyl-3-aminopropionate is C{amine['C']}H{amine['H']}N{amine['N']}O{amine['O']}.")

    # Step 3: The reaction is S+ + amine -> B+ + H2O.
    # The byproduct is water (H2O).
    water = {'C': 0, 'H': 2, 'O': 1, 'N': 0}
    print(f"Step 3: The reaction produces the product cation B+ and a water molecule (H2O).")

    # Step 4: Calculate the atomic composition of the product cation B+.
    # B+ = S+ + amine - H2O
    b_cation_c = s_cation['C'] + amine['C'] - water['C']
    b_cation_h = s_cation['H'] + amine['H'] - water['H']
    b_cation_n = s_cation['N'] + amine['N'] - water['N']
    b_cation_o = s_cation['O'] + amine['O'] - water['O']

    print("\nCalculating the formula for the cation of compound B:")
    print(f"Carbon atoms: {s_cation['C']} (from S+) + {amine['C']} (from amine) - {water['C']} (from H2O) = {b_cation_c}")
    print(f"Hydrogen atoms: {s_cation['H']} (from S+) + {amine['H']} (from amine) - {water['H']} (from H2O) = {b_cation_h}")
    print(f"Nitrogen atoms: {s_cation['N']} (from S+) + {amine['N']} (from amine) - {water['N']} (from H2O) = {b_cation_n}")
    print(f"Oxygen atoms: {s_cation['O']} (from S+) + {amine['O']} (from amine) - {water['O']} (from H2O) = {b_cation_o}")

    # Output the final molecular formula for the organic cation of B.
    final_formula = f"C{b_cation_c}H{b_cation_h}N{b_cation_n}O{b_cation_o}"
    print(f"\nThe molecular formula of the organic cation of compound B is {final_formula}.")

calculate_molecular_formula_B()