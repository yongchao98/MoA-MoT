import operator

def find_optimal_material():
    """
    This function compares materials based on their radiation length (X0)
    to find the optimal choice for a particle detector cooling system.
    A longer radiation length is better as it minimizes particle interaction.
    """
    # Data for materials: Radiation Length (X0) in cm.
    # Source: Particle Data Group (PDG)
    materials = {
        "Titanium": 3.56,
        "Aluminium": 8.9,
        "316 Stainless Steel": 1.76, # Approx. value for iron/steel
        "Copper": 1.43,
        "Nickel": 1.42
    }

    print("Comparing materials based on Radiation Length (X0).")
    print("A higher value means less particle interaction, which is better for a detector.\n")

    # Find the material with the maximum radiation length
    optimal_material = max(materials.items(), key=operator.itemgetter(1))

    # Print the comparison
    for material, x0 in materials.items():
        print(f"- {material}: Radiation Length = {x0} cm")

    print(f"\nConclusion:")
    print(f"For maximizing the unique parameter (Radiation Length), {optimal_material[0]} is the optimum choice.")
    print(f"Final Equation: Maximum Radiation Length = {optimal_material[1]} cm")


if __name__ == "__main__":
    find_optimal_material()