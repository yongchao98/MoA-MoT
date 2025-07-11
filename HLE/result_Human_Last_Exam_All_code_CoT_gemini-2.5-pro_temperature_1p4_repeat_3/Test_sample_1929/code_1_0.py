def calculate_molar_mass():
    """
    Identifies the product of butadiene and 1,1-dichloro-2,2-difluoroethene,
    and calculates its molar mass.
    """
    # The reaction is a Diels-Alder [4+2] cycloaddition.
    # Reactants:
    # 1,3-butadiene (C4H6)
    # 1,1-dichloro-2,2-difluoroethene (C2Cl2F2)
    # The reaction forms a six-membered ring.
    # Product: 4,4-dichloro-5,5-difluorocyclohex-1-ene
    # Chemical Formula: C6H6Cl2F2

    product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"
    product_formula = "C6H6Cl2F2"

    print(f"The product of the reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is {product_name}.")
    print(f"The chemical formula of the product is {product_formula}.")
    print("\nThe molar mass of the product is calculated below.")

    # Standard atomic weights (g/mol)
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'Cl': 35.453,
        'F': 18.998
    }

    # Atom counts in the molecule C6H6Cl2F2
    molecule_counts = {
        'C': 6,
        'H': 6,
        'Cl': 2,
        'F': 2
    }

    total_mass = 0
    equation_terms = []
    
    # Build the equation string parts and calculate total mass
    # The order is just for display purposes
    for element in ['C', 'H', 'Cl', 'F']:
        count = molecule_counts[element]
        weight = atomic_weights[element]
        contribution = count * weight
        total_mass += contribution
        # Create a string for each part of the equation, e.g., "(6 * 12.011)"
        equation_terms.append(f"({count} * {weight})")

    # Print the full equation showing each number
    equation_str = " + ".join(equation_terms)
    print(f"Molar Mass = {equation_str}")
    
    # Calculate and print the value of each term
    term_values = []
    for element in ['C', 'H', 'Cl', 'F']:
        count = molecule_counts[element]
        weight = atomic_weights[element]
        term_values.append(f"{(count * weight):.3f}")

    intermediate_str = " + ".join(term_values)
    print(f"Molar Mass = {intermediate_str}")
    
    # Print the final result
    print(f"Molar Mass = {total_mass:.3f} g/mol")

# Execute the function
calculate_molar_mass()