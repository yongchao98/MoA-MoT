def find_optimal_material():
    """
    Determines the optimal material for a particle detector's cooling system
    by maximizing the radiation length.
    """

    print("Analyzing materials for a particle detector's cooling system.")
    print("The unique and most critical parameter for this application is minimizing interaction with the particles being detected.")
    print("This is achieved by selecting a material with the longest possible 'radiation length' (X0).")
    print("A longer radiation length means the material is more 'transparent' to particles, causing less scattering and measurement error.")
    print("-" * 30)

    # Data: Material candidates and their approximate radiation lengths in cm.
    # A longer radiation length is better for a particle detector.
    materials = {
        "A. Titanium": {"radiation_length_cm": 3.6},
        "B. Aluminium": {"radiation_length_cm": 8.9},
        "C. 316 Stainless Steel": {"radiation_length_cm": 1.75}, # Using value for Iron/Steel
        "D. Copper": {"radiation_length_cm": 1.43},
        "E. Nickle": {"radiation_length_cm": 1.42}
    }

    # Find the material with the maximum radiation length
    best_material_name = None
    max_radiation_length = -1

    print("Comparing the radiation lengths (X0) of the candidate materials:")
    for name, properties in materials.items():
        length = properties["radiation_length_cm"]
        print(f"{name}: X0 = {length} cm")
        if length > max_radiation_length:
            max_radiation_length = length
            best_material_name = name

    print("-" * 30)
    print("Conclusion:")
    print(f"The optimum choice is the material with the maximum radiation length, which is {best_material_name}.")
    print(f"For pipes of similar dimensions, Aluminium will cause the least disturbance to the particle measurements.")

find_optimal_material()
<<<B>>>