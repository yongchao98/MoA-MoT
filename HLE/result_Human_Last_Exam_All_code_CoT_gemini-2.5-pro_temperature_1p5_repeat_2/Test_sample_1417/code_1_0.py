import operator

def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling tubes
    by finding the material with the longest radiation length.
    """
    # Data for materials: Radiation Length in cm.
    # Radiation Length (X0) is the key parameter. A longer X0 means less
    # particle interaction, which is crucial for a particle detector.
    # Formula: X0 (cm) = X0 (g/cm^2) / density (g/cm^3)
    materials = {
        'A. Titanium': 3.56,
        'B. Aluminium': 8.9,
        'C. 316 Stainless Steel': 1.74, # Approximate, similar to Iron
        'D. Copper': 1.43,
        'E. Nickle': 1.42
    }

    print("Comparing materials based on the unique parameter for a particle detector: Radiation Length (X0).")
    print("A longer radiation length is better as it minimizes interference with the particles being measured.")
    print("-" * 50)

    # Print the radiation length for each material
    for material, x0 in materials.items():
        print(f"Material: {material}, Radiation Length (X0) = {x0} cm")

    # Find the material with the maximum radiation length
    best_material_key = max(materials.items(), key=operator.itemgetter(1))[0]
    
    print("-" * 50)
    print(f"Conclusion: {best_material_key} has the longest radiation length, making it the optimum choice.")
    print("This minimizes the material's interaction with particles, which is the unique and critical requirement for this application.")

find_optimal_material()