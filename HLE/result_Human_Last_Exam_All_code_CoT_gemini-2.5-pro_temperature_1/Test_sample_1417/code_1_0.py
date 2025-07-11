def find_optimal_material():
    """
    This function compares different materials based on their radiation length
    to find the optimal choice for a particle detector cooling system.
    """
    materials = {
        'Titanium': {'Z': 22, 'X0_cm': 3.56},
        'Aluminium': {'Z': 13, 'X0_cm': 8.9},
        '316 Stainless Steel': {'Z': '~26 (avg)', 'X0_cm': 1.74}, # Approximation based on Iron
        'Copper': {'Z': 29, 'X0_cm': 1.43},
        'Nickle': {'Z': 28, 'X0_cm': 1.42}
    }

    print("Comparing materials based on Radiation Length (X0).")
    print("A longer radiation length is better as it minimizes particle interaction.")
    print("-" * 60)
    print(f"{'Material':<25} | {'Atomic Number (Z)':<20} | {'Radiation Length (cm)':<20}")
    print("-" * 60)

    optimal_material = None
    max_x0 = -1

    for name, properties in materials.items():
        z_val = properties['Z']
        x0_val = properties['X0_cm']
        print(f"{name:<25} | {str(z_val):<20} | {x0_val:<20.2f}")
        if x0_val > max_x0:
            max_x0 = x0_val
            optimal_material = name

    print("-" * 60)
    print(f"\nThe material with the longest radiation length is {optimal_material}.")
    print("Therefore, it is the optimum choice to minimize the impact on particle measurements.")

find_optimal_material()