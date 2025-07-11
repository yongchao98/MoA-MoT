def calculate_molecular_formula():
    """
    Calculates the molecular formula of the starting material,
    (3,4-dihydro-2H-pyrrol-5-yl)proline, by summing the atoms of its constituent parts.
    """
    # The starting material is N-(3,4-dihydro-2H-pyrrol-5-yl)pyrrolidine-2-carboxylic acid.
    # We break it down into two parts for atom counting.

    # Part 1: The N-substituted proline core (pyrrolidine-2-carboxylic acid part).
    # - Carbons: 4 in the pyrrolidine ring + 1 in the -COOH group = 5
    # - Hydrogens: 1 on C2, 2 on C3, 2 on C4, 2 on C5, and 1 on the -COOH group = 8
    # - Nitrogens: 1 in the ring.
    # - Oxygens: 2 in the -COOH group.
    proline_part = {'C': 5, 'H': 8, 'N': 1, 'O': 2}

    # Part 2: The (3,4-dihydro-2H-pyrrol-5-yl) substituent.
    # This is a cyclic imine attached via its C5' carbon.
    # - Carbons: 4 in the ring.
    # - Hydrogens: 2 on C2', 2 on C3', 2 on C4' = 6
    # - Nitrogens: 1 in the ring.
    # - Oxygens: 0.
    substituent_part = {'C': 4, 'H': 6, 'N': 1, 'O': 0}

    # Sum the atoms from both parts
    total_C = proline_part['C'] + substituent_part['C']
    total_H = proline_part['H'] + substituent_part['H']
    total_N = proline_part['N'] + substituent_part['N']
    total_O = proline_part['O'] + substituent_part['O']

    print("Calculating the molecular formula of the starting material:")
    print("\nStep 1: Analyze the N-substituted proline core.")
    print(f"Formula: C{proline_part['C']}H{proline_part['H']}N{proline_part['N']}O{proline_part['O']}")

    print("\nStep 2: Analyze the (3,4-dihydro-2H-pyrrol-5-yl) substituent.")
    print(f"Formula: C{substituent_part['C']}H{substituent_part['H']}N{substituent_part['N']}")

    print("\nStep 3: Sum the atoms from both parts.")
    print(f"Total Carbons (C)   = {proline_part['C']} + {substituent_part['C']} = {total_C}")
    print(f"Total Hydrogens (H) = {proline_part['H']} + {substituent_part['H']} = {total_H}")
    print(f"Total Nitrogens (N) = {proline_part['N']} + {substituent_part['N']} = {total_N}")
    print(f"Total Oxygens (O)   = {proline_part['O']} + {substituent_part['O']} = {total_O}")

    print(f"\nThe final molecular formula of the starting material is C{total_C}H{total_H}N{total_N}O{total_O}.")

calculate_molecular_formula()
<<<C9H14N2O2>>>