def calculate_compound_b_formula():
    """
    This script calculates the molecular formula of compound B based on the reaction scheme.
    """

    # Step 1: Define the atomic compositions of the species involved in the formation of Cation B.
    # The precursor xanthenium cation, as deduced from the structure of A.
    precursor_cation = {'C': 23, 'H': 21, 'O': 5}
    # The amine reactant, methyl-3-aminopropionate (H2N-CH2CH2COOCH3).
    amine = {'C': 4, 'H': 9, 'N': 1, 'O': 2}
    # The byproduct of the condensation reaction.
    water = {'H': 2, 'O': 1}

    # Print the reaction to form the cation of B
    print("The chemical equation for the formation of the cation of compound B is:")
    # Reactants
    precursor_str = f"C{precursor_cation['C']}H{precursor_cation['H']}O{precursor_cation['O']}(+)"
    amine_str = f"C{amine['C']}H{amine['H']}N{amine['N']}O{amine['O']}"
    # Products
    water_str = f"H{water['H']}O{water['O']}"
    
    # Step 2: Calculate the atomic composition of the cation of B.
    # Cation B = Precursor + Amine - Water
    cation_B = {}
    all_reactants = precursor_cation.copy()
    for element, count in amine.items():
        all_reactants[element] = all_reactants.get(element, 0) + count

    cation_B = all_reactants.copy()
    for element, count in water.items():
        cation_B[element] -= count
        if cation_B[element] == 0:
            del cation_B[element]

    cation_B_str = f"C{cation_B['C']}H{cation_B['H']}N{cation_B['N']}O{cation_B['O']}(+)"
    print(f"{precursor_str} + {amine_str} -> {cation_B_str} + {water_str}")
    print("-" * 20)

    # Step 3: Add the counter-ion (BF4-) to get the final formula of compound B.
    counter_ion = {'B': 1, 'F': 4}
    
    # Print the salt formation step
    print("The cation then forms a salt with the tetrafluoroborate anion:")
    counter_ion_str = f"B{counter_ion['B']}F{counter_ion['F']}(-)"
    
    final_compound = cation_B.copy()
    for element, count in counter_ion.items():
        final_compound[element] = final_compound.get(element, 0) + count

    # Step 4: Format the final molecular formula string in the standard order (C, H, B, F, N, O).
    final_formula_str = (f"C{final_compound['C']}"
                         f"H{final_compound['H']}"
                         f"B"  # Count is 1
                         f"F{final_compound['F']}"
                         f"N"  # Count is 1
                         f"O{final_compound['O']}")
    
    print(f"{cation_B_str} + {counter_ion_str} -> {final_formula_str}")
    print("-" * 20)
    print(f"The molecular formula of compound B is: {final_formula_str}")

calculate_compound_b_formula()