import operator

def find_optimal_material():
    """
    Analyzes materials to find the one with the longest radiation length,
    which is the optimal choice for minimizing particle interaction in a detector.
    """
    # Data for materials: Radiation Length (X0) in cm.
    # A longer radiation length is better for a particle detector.
    materials = {
        'Titanium': 3.56,
        'Aluminium': 8.9,
        '316 Stainless Steel': 1.76,
        'Copper': 1.43,
        'Nickle': 1.42
    }

    print("The unique parameter for a particle detector cooling system is minimizing particle interaction.")
    print("This is quantified by the material's radiation length (X0). A longer radiation length is better.")
    print("\nComparing the radiation lengths (in cm) of the candidate materials:")

    for material, x0 in materials.items():
        print(f"- {material}: {x0} cm")

    # Find the material with the maximum radiation length
    # The `operator.itemgetter(1)` function gets the value (the second item) from each dictionary item
    optimal_material, max_x0 = max(materials.items(), key=operator.itemgetter(1))

    print(f"\nThe material with the longest radiation length is {optimal_material} with X0 = {max_x0} cm.")
    print(f"Therefore, {optimal_material} is the optimum choice.")

if __name__ == "__main__":
    find_optimal_material()