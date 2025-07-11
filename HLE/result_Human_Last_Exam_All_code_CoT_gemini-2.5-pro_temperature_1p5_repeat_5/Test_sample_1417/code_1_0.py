def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    by finding the material with the maximum radiation length.
    """
    # The unique parameter for a material in a particle detector is its interaction
    # with the particles being measured. This is minimized by maximizing the
    # material's 'Radiation Length' (X0). A longer X0 is better.

    # Data from the Particle Data Group (PDG).
    # Format: { 'Material Name': [Atomic Number (Z), Radiation Length (X0) in cm] }
    # 316 Stainless steel is approximated by its main component, Iron (Fe).
    materials_data = {
        'A. Titanium': [22, 3.56],
        'B. Aluminium': [13, 8.897],
        'C. 316 Stainless Steel': [26, 1.757],
        'D. Copper': [29, 1.436],
        'E. Nickle': [28, 1.424]
    }

    # Initialize variables to find the best material
    optimal_material_name = None
    max_radiation_length = -1.0

    print("Goal: Find the material with the maximum Radiation Length (X0).")
    print("A longer radiation length means less interaction with particles, which is ideal for a particle detector.\n")
    print("Comparing Materials:")
    print(f"{'Choice':<25} {'Atomic Number (Z)':<20} {'Radiation Length (X0/cm)':<25}")
    print("-" * 70)

    # Iterate through the materials to find the one with the maximum X0
    for name, properties in materials_data.items():
        z_number = properties[0]
        radiation_length = properties[1]
        
        print(f"{name:<25} {z_number:<20} {radiation_length:<25.3f}")

        if radiation_length > max_radiation_length:
            max_radiation_length = radiation_length
            optimal_material_name = name
    
    print("-" * 70)
    print(f"\nConclusion: The material with the longest radiation length is '{optimal_material_name}'.")
    print("Therefore, it is the optimum choice for minimizing particle interaction in the detector.")


if __name__ == "__main__":
    find_optimal_material()