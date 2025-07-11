def calculate_product_formula():
    """
    Calculates the molecular formula of the final product from the given reaction scheme.
    """
    # Step 1: Calculate the formula of the intermediate
    # Reactant 1: 2-aminothiazole
    reactant1 = {'C': 3, 'H': 4, 'N': 2, 'S': 1, 'Cl': 0, 'O': 0}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate
    reactant2 = {'C': 6, 'H': 9, 'N': 0, 'S': 0, 'Cl': 1, 'O': 3}

    # Byproducts from the first reaction (condensation)
    hcl = {'C': 0, 'H': 1, 'Cl': 1, 'N': 0, 'O': 0, 'S': 0}
    h2o = {'C': 0, 'H': 2, 'O': 1, 'N': 0, 'Cl': 0, 'S': 0}

    # Calculate the intermediate's formula
    intermediate = {}
    all_elements = set(reactant1.keys()) | set(reactant2.keys())
    for element in all_elements:
        intermediate[element] = reactant1.get(element, 0) + reactant2.get(element, 0) \
                                - hcl.get(element, 0) - h2o.get(element, 0)

    # Intermediate formula is C9H10N2O2S

    # Step 2: Calculate the formula of the final product
    # The reaction is an amidation of the ester.
    # Group removed from intermediate: ethoxy group (-OEt)
    ethoxy_group = {'C': 2, 'H': 5, 'O': 1}
    # Group added: benzylamino group (-NH-CH2-Ph)
    benzylamino_group = {'C': 7, 'H': 8, 'N': 1}

    # Calculate the final product's formula
    final_product = intermediate.copy()
    for element, count in ethoxy_group.items():
        final_product[element] -= count
    for element, count in benzylamino_group.items():
        final_product[element] = final_product.get(element, 0) + count

    # Print the counts for each element in the final product
    print("The molecular formula of the final product is determined by the count of each element:")
    print(f"Carbon (C): {final_product['C']}")
    print(f"Hydrogen (H): {final_product['H']}")
    print(f"Nitrogen (N): {final_product['N']}")
    print(f"Oxygen (O): {final_product['O']}")
    print(f"Sulfur (S): {final_product['S']}")
    
    # Format the final molecular formula string
    formula = (f"C{final_product['C']}"
               f"H{final_product['H']}"
               f"N{final_product['N']}"
               f"O{final_product['O'] if final_product['O'] > 1 else ''}"
               f"S{final_product['S'] if final_product['S'] > 1 else ''}")
    
    print(f"\nThe final molecular formula is: {formula}")

calculate_product_formula()