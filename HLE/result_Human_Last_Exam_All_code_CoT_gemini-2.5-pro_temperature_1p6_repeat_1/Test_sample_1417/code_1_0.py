def find_optimal_material():
    """
    This script identifies the optimal material for a particle detector's cooling system
    based on the unique requirement of maximizing particle transparency.

    This transparency is quantified by the material's radiation length (X0). A longer
    radiation length is better as it minimizes the interaction and scattering of particles
    that the detector is trying to measure.
    """

    # Radiation lengths (X0) in cm for various materials. A higher value is better.
    # Source: Particle Data Group (PDG)
    materials = {
        "Titanium": 3.56,
        "Aluminium": 8.9,
        "316 Stainless Steel": 1.76, # Approximated using Iron, its main component
        "Copper": 1.43,
        "Nickle": 1.47
    }

    # Find the material with the maximum radiation length
    optimal_material = max(materials, key=materials.get)
    max_x0 = materials[optimal_material]

    print("Analyzing materials based on the unique requirement for a particle detector:")
    print("The goal is to maximize the radiation length (X0) to minimize particle interaction.\n")

    print("Material Radiation Lengths (X0 in cm):")
    print("---------------------------------------")
    for material, x0 in materials.items():
        print(f"- {material}: {x0} cm")
    print("---------------------------------------\n")

    print(f"Conclusion: To maximize the radiation length, the optimum choice is {optimal_material}.")
    print(f"Final Equation: {optimal_material} has the max radiation length of {max_x0} cm.")

if __name__ == '__main__':
    find_optimal_material()
<<<B>>>