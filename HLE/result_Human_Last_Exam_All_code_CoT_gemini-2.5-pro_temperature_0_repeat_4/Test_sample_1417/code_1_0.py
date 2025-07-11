def find_optimal_material():
    """
    Analyzes materials to find the best choice for a particle detector's
    cooling system based on minimizing particle interaction.
    """
    materials = {
        'A. Titanium': {'Z': 22, 'density_g_cm3': 4.51},
        'B. Aluminium': {'Z': 13, 'density_g_cm3': 2.70},
        'C. 316 Stainless Steel': {'Z': '~26 (Fe)', 'density_g_cm3': 8.00},
        'D. Copper': {'Z': 29, 'density_g_cm3': 8.96},
        'E. Nickle': {'Z': 28, 'density_g_cm3': 8.90}
    }

    print("Comparing materials based on properties relevant to particle interaction.")
    print("The goal is to minimize the material budget by choosing a material with low Atomic Number (Z) and low density.")
    print("-" * 70)
    print(f"{'Material':<25} | {'Approx. Atomic Number (Z)':<30} | {'Density (g/cm^3)'}")
    print("-" * 70)

    for name, properties in materials.items():
        z_val = properties['Z']
        density_val = properties['density_g_cm3']
        print(f"{name:<25} | {str(z_val):<30} | {density_val:.2f}")

    print("-" * 70)
    print("\nConclusion:")
    print("Aluminium has the lowest Atomic Number and the lowest density of the choices.")
    print("This means it will interact the least with particles passing through it, making it the optimal choice for preserving the accuracy of particle measurements.")

find_optimal_material()