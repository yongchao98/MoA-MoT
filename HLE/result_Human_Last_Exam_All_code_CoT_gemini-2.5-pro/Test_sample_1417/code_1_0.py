def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling system
    based on the unique parameter of minimizing particle interaction.

    In a particle detector, the primary goal is to measure particles originating
    from a collision with minimal interference. Materials within the detector's
    path can scatter or stop these particles, corrupting the measurement. This
    effect is primarily dependent on the material's atomic number (Z). A lower
    atomic number leads to less interaction.

    This script compares the candidate materials based on their atomic number
    to find the one that would be most "transparent" to particles.
    """

    # Data for candidate materials: Atomic Number (Z)
    # For alloys like Stainless Steel, we use the Z of the primary element (Iron).
    materials = {
        "Titanium": {"Z": 22, "Choice": "A"},
        "Aluminium": {"Z": 13, "Choice": "B"},
        "316 Stainless Steel": {"Z": 26, "Choice": "C"},
        "Copper": {"Z": 29, "Choice": "D"},
        "Nickle": {"Z": 28, "Choice": "E"}
    }

    # Initialize variables to find the material with the minimum Z
    optimal_material = None
    min_z = float('inf')

    print("Evaluating materials based on minimizing particle interaction (lowest Atomic Number Z):")
    for name, properties in materials.items():
        z_value = properties["Z"]
        print(f"- {name} (Z = {z_value})")
        if z_value < min_z:
            min_z = z_value
            optimal_material = name

    print("\nConclusion:")
    print(f"The unique parameter for a particle detector is minimizing particle interaction.")
    print(f"This is achieved by selecting the material with the lowest atomic number (Z).")
    print(f"Based on the data, {optimal_material} has the lowest atomic number (Z = {min_z}).")
    print(f"Therefore, {optimal_material} is the optimum choice.")

find_optimal_material()