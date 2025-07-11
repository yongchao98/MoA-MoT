def calculate_molecular_weight():
    """
    Calculates and prints the molecular weight of Compound A, which is
    (phenylamino)(3-hydroxypyridin-2-yl)acetonitrile.
    """
    # Standard atomic weights (g/mol)
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Molecular formula of Compound A
    formula = {'C': 13, 'H': 11, 'N': 3, 'O': 1}

    molecular_weight = 0
    calculation_parts = []

    # Iterate through the elements in the formula to build the calculation string and sum the weights
    for element, count in formula.items():
        weight = atomic_weights[element]
        molecular_weight += count * weight
        calculation_parts.append(f"{count} * {weight} ({element})")

    # Join the parts into a final equation string
    calculation_equation = " + ".join(calculation_parts)

    print("The final product, Compound A, is (phenylamino)(3-hydroxypyridin-2-yl)acetonitrile.")
    print("Its molecular formula is C13H11N3O.")
    print("\nThe molecular weight is calculated as follows:")
    print(f"{calculation_equation} = {molecular_weight:.3f} g/mol")

# Execute the function to display the result
calculate_molecular_weight()