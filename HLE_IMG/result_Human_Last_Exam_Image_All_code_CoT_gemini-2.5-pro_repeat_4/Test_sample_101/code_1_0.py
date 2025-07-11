def calculate_molecular_weight():
    """
    Calculates and prints the molecular weight of Compound A.
    The reaction is an imine formation followed by a Strecker-type reaction.
    1. 3-hydroxy-pyridine-2-carbaldehyde + Aniline -> Imine intermediate
    2. Imine intermediate + NaCN -> Compound A (an alpha-aminonitrile)
    Compound A is 2-(anilino(cyano)methyl)pyridin-3-ol.
    Its molecular formula is C13H11N3O.
    """

    # Standard atomic weights (g/mol)
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Atom counts for Compound A (C13H11N3O)
    atom_counts = {
        'C': 13,
        'H': 11,
        'N': 3,
        'O': 1
    }

    # Calculate the molecular weight
    molecular_weight = 0
    calculation_parts = []
    
    # Iterate through the elements to build the calculation string and sum the weights
    for element in ['C', 'H', 'N', 'O']:
        count = atom_counts[element]
        weight = atomic_weights[element]
        term_weight = count * weight
        molecular_weight += term_weight
        calculation_parts.append(f"({count} * {weight})")

    calculation_equation = " + ".join(calculation_parts)

    print("The final product, Compound A, is 2-(anilino(cyano)methyl)pyridin-3-ol.")
    print("The molecular formula is C13H11N3O.")
    print("\nMolecular Weight Calculation:")
    print(f"MW = C*13 + H*11 + N*3 + O*1")
    print(f"MW = {calculation_equation}")
    print(f"MW = {molecular_weight:.3f} g/mol")

calculate_molecular_weight()