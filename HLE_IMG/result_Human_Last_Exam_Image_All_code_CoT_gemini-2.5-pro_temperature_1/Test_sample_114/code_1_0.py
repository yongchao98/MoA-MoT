def calculate_exact_mass():
    """
    Calculates the exact mass of chemical compounds based on their molecular formulas
    to verify the High-Resolution Mass Spectrometry (HRMS) data provided.
    """
    # Exact masses of the most abundant isotopes
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915
    }

    # Molecular formulas derived from the problem description's HRMS data
    formulas = {
        "Product A": {'C': 9, 'H': 13, 'N': 1, 'O': 2},
        "Product B": {'C': 10, 'H': 13, 'N': 1, 'O': 2}
    }

    print("Verifying molecular formulas by calculating exact mass:\n")

    for product, elements in formulas.items():
        total_mass = 0
        equation_parts = []
        # Build the equation string part by part
        for element, count in elements.items():
            mass = atomic_masses[element]
            total_mass += count * mass
            equation_parts.append(f"{count} * {mass:.6f}") # Representing the calculation for each element
        
        formula_str = ''.join([f'{k}{v}' for k, v in elements.items()])
        print(f"Calculation for {product} ({formula_str}):")
        
        # Print the full equation showing each term
        equation_str_full = f"Exact Mass = ({equation_parts[0]})_C + ({equation_parts[1]})_H + ({equation_parts[2]})_N + ({equation_parts[3]})_O"
        # Since the order is fixed from the dictionary, we can use indices
        print(equation_str_full)
        
        # Print the final result
        print(f"Resulting Exact Mass for {product}: {total_mass:.5f}\n")

# Execute the function to show the calculations
calculate_exact_mass()