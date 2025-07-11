def solve_material_choice():
    """
    Analyzes materials for a particle detector cooling system.

    The unique parameter for a particle detector application, compared to a spacecraft,
    is the need to minimize the interaction between the material and the particles being
    measured. This is quantified by the material's radiation length (X0). A longer
    radiation length is better as it signifies less interaction.

    This script compares the radiation lengths of the candidate materials to find the
    optimum choice.
    """
    materials = {
        'A. Titanium': {'symbol': 'Ti', 'Z': 22, 'radiation_length_cm': 3.56},
        'B. Aluminium': {'symbol': 'Al', 'Z': 13, 'radiation_length_cm': 8.9},
        'C. 316 Stainless Steel': {'symbol': 'Fe (approx)', 'Z': 26, 'radiation_length_cm': 1.76},
        'D. Copper': {'symbol': 'Cu', 'Z': 29, 'radiation_length_cm': 1.43},
        'E. Nickel': {'symbol': 'Ni', 'Z': 28, 'radiation_length_cm': 1.42},
    }

    print("--- Material Analysis for Particle Detector Cooling System ---")
    print("The unique and critical parameter to maximize is the Radiation Length (X0).")
    print("A longer radiation length means less interference with particle measurements.\n")
    print(f"{'Choice':<25} | {'Atomic Number (Z)':<20} | {'Radiation Length (cm)':<20}")
    print("-" * 75)

    best_material = None
    max_radiation_length = -1

    for name, properties in materials.items():
        print(f"{name:<25} | {properties['Z']:<20} | {properties['radiation_length_cm']:<20.2f}")
        if properties['radiation_length_cm'] > max_radiation_length:
            max_radiation_length = properties['radiation_length_cm']
            best_material = name

    print("-" * 75)
    print(f"\nConclusion: {best_material} has the longest radiation length ({max_radiation_length:.2f} cm).")
    print("Therefore, Aluminium is the optimum choice to minimize the material's impact on the particle measurements.")

solve_material_choice()
<<<B>>>