def solve_molecular_formula():
    """
    This script calculates the molecular formula of compound B based on the provided reaction scheme.
    """
    # Step 1: Define the atomic composition of the starting xanthylium cation.
    # From structural analysis of the image (9-phenyl-tetramethoxyxanthylium cation),
    # its formula is C23H21O5.
    start_cation = {'C': 23, 'H': 21, 'O': 5}

    # Step 2: Define the atomic composition of the reactant, methyl-3-aminopropionate (H2N-CH2-CH2-COOCH3).
    # Its formula is C4H9NO2.
    amine = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

    # Step 3: Define the atomic composition of the byproduct, water (H2O), which is eliminated.
    water = {'H': 2, 'O': 1}

    # Step 4: Calculate the composition of the product cation B+.
    # The reaction is: Start_Cation + Amine -> Product_Cation + Water
    product_cation = {}
    all_elements_in_reaction = set(start_cation.keys()) | set(amine.keys()) | set(water.keys())
    for element in sorted(list(all_elements_in_reaction)):
        count = start_cation.get(element, 0) + amine.get(element, 0) - water.get(element, 0)
        if count > 0:
            product_cation[element] = count
            
    # Step 5: Define the counter-ion, tetrafluoroborate (BF4-), added in the final step.
    anion = {'B': 1, 'F': 4}

    # Step 6: Calculate the final molecular formula of compound B (the salt).
    final_compound = product_cation.copy()
    for element, count in anion.items():
        final_compound[element] = final_compound.get(element, 0) + count

    # Step 7: Format the output as requested.
    # The prompt asks to output each number in the final equation.
    # We will show the elemental composition of the final neutral compound.
    print("The final compound B is a salt composed of the organic cation and the BF4- anion.")
    print("The elemental composition is calculated as follows:")
    
    # Generate the equation string showing the number for each atom
    equation_parts = []
    # Standard order for formula: C, H, then alphabetical
    element_order = ['C', 'H'] + sorted([el for el in final_compound.keys() if el not in ['C', 'H']])
    
    for el in element_order:
        if el in final_compound:
            equation_parts.append(f"{el}{final_compound[el]}")

    print(" + ".join(equation_parts))

    # Generate the final molecular formula string in Hill system format.
    final_formula_str = ""
    if 'C' in final_compound:
        final_formula_str += f"C{final_compound['C']}"
    if 'H' in final_compound:
        final_formula_str += f"H{final_compound['H']}"
    
    # Add other elements alphabetically
    other_elements = sorted([el for el in final_compound if el not in ['C', 'H']])
    for el in other_elements:
        count = final_compound[el]
        final_formula_str += f"{el}{count if count > 1 else ''}"

    print(f"\nThe molecular formula of compound B is: {final_formula_str}")

solve_molecular_formula()