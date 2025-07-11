def get_diels_alder_product():
    """
    Determines the product of the Diels-Alder reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the reaction equation.
    """

    # Helper function to format a chemical formula from a dictionary of atoms
    def format_formula(formula_dict):
        # Standard order for organic compounds: C, H, then others alphabetically
        elements = ['C', 'H'] + sorted([el for el in formula_dict if el not in ['C', 'H']])
        formula_str = ""
        for el in elements:
            if el in formula_dict:
                count = formula_dict[el]
                formula_str += el
                if count > 1:
                    formula_str += str(count)
        return formula_str

    # 1. Define the reactants with their names and chemical formulas
    diene = {
        'name': 'Buta-1,3-diene',
        'formula': {'C': 4, 'H': 6}
    }

    dienophile = {
        'name': '1,1-dichloro-2,2-difluoroethene',
        'formula': {'C': 2, 'Cl': 2, 'F': 2}
    }

    # 2. Determine the product by adding the atoms of the reactants
    product_formula = {}
    all_elements = set(diene['formula'].keys()) | set(dienophile['formula'].keys())
    for element in all_elements:
        product_formula[element] = diene['formula'].get(element, 0) + dienophile['formula'].get(element, 0)

    # 3. Define the product name based on IUPAC nomenclature
    product = {
        'name': '4,4-dichloro-5,5-difluorocyclohex-1-ene',
        'formula': product_formula
    }

    # 4. Format all chemical formulas into strings for printing
    diene_formula_str = format_formula(diene['formula'])
    dienophile_formula_str = format_formula(dienophile['formula'])
    product_formula_str = format_formula(product['formula'])

    # 5. Print the details and the final balanced equation
    print("The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder cycloaddition.")
    print("\n--- Reactants ---")
    print(f"Diene: {diene['name']} ({diene_formula_str})")
    print(f"Dienophile: {dienophile['name']} ({dienophile_formula_str})")
    
    print("\n--- Product ---")
    print(f"Product: {product['name']} ({product_formula_str})")

    print("\n--- Final Equation ---")
    # The stoichiometric coefficients are all 1
    # The output string contains all numbers from the chemical formulas
    print(f"1 {diene_formula_str} + 1 {dienophile_formula_str} -> 1 {product_formula_str}")

# Execute the function
get_diels_alder_product()