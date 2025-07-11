def calculate_molecular_weight():
    """
    This function identifies Compound A, calculates its molecular formula,
    and computes its molecular weight, showing each step of the calculation.
    """
    # Define atomic masses for the elements involved
    atomic_mass = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999}

    # Compound A is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.
    # Let's determine its molecular formula by counting the atoms:
    # Pyridine ring: 5C, 3H, 1N
    # OH group: 1O, 1H
    # Phenylamino group: 6C, 5H (on ring) + 1H (on N), 1N
    # Acetonitrile part (excluding the N already counted in phenylamino): 1C (alpha), 1H (alpha), 1C (cyano), 1N (cyano)
    # Total C = 5 + 6 + 1 + 1 = 13
    # Total H = 3 + 1 + 5 + 1 + 1 = 11
    # Total N = 1 + 1 + 1 = 3
    # Total O = 1
    # Formula: C13H11N3O
    
    num_atoms = {'C': 13, 'H': 11, 'N': 3, 'O': 1}

    print("The final product, Compound A, is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.")
    print(f"Molecular Formula: C{num_atoms['C']}H{num_atoms['H']}N{num_atoms['N']}O{num_atoms['O']}")
    print("\nCalculation of the molecular weight:")

    total_mw = 0
    calculation_parts = []
    
    # Calculate and display the mass contribution of each element
    for element, count in num_atoms.items():
        mass = count * atomic_mass[element]
        total_mw += mass
        print(f"Contribution from {element}: {count} * {atomic_mass[element]} = {mass:.3f}")
        calculation_parts.append(f"({count} * {atomic_mass[element]})")

    # Display the final summation equation
    print("\nThe final calculation is:")
    final_equation = " + ".join(calculation_parts)
    print(f"{final_equation} = {total_mw:.3f} g/mol")
    
    print(f"\nThe molecular weight of Compound A is approximately {total_mw:.3f} g/mol.")

calculate_molecular_weight()