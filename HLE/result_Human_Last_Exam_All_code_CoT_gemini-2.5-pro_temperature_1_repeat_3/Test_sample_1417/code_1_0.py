def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    based on the unique parameter of minimizing particle interaction.
    """
    materials = [
        {'name': 'Aluminium', 'Z': 13, 'density_g_cm3': 2.70},
        {'name': 'Titanium', 'Z': 22, 'density_g_cm3': 4.51},
        {'name': 'Nickel', 'Z': 28, 'density_g_cm3': 8.91},
        {'name': 'Copper', 'Z': 29, 'density_g_cm3': 8.96},
        {'name': '316 Stainless Steel', 'Z': 26, 'density_g_cm3': 8.00}, # Using Iron (Fe) as the primary component
    ]

    print("Evaluating materials for a particle detector cooling system.")
    print("The unique and most critical parameter is minimizing particle interaction.")
    print("This is achieved by selecting a material with a low atomic number (Z) and low density (rho).\n")

    # Sort materials by Z number to find the minimum
    # In this case, sorting by Z is sufficient as the material with the lowest Z also has the lowest density.
    optimal_material = min(materials, key=lambda x: (x['Z'], x['density_g_cm3']))

    print("Material properties:")
    for mat in materials:
        print(f"- {mat['name']:<22} Atomic Number (Z): {mat['Z']:<4} Density (g/cm^3): {mat['density_g_cm3']}")

    print("\n--- Conclusion ---")
    print(f"Comparing the candidates, {optimal_material['name']} has the lowest atomic number (Z={optimal_material['Z']}) and the lowest density ({optimal_material['density_g_cm3']} g/cm^3).")
    print(f"This results in the longest radiation length, causing the least interference with particle measurements.")
    print(f"Therefore, {optimal_material['name']} is the optimum choice.")

find_optimal_material()