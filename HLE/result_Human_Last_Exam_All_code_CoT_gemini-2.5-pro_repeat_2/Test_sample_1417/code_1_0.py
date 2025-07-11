import sys

def solve_detector_material_choice():
    """
    This script determines the optimum material for a particle detector's cooling tubes
    by finding the material with the longest radiation length (X₀), which minimizes
    interaction with the particles being measured.
    """

    # Data for materials: Radiation Length (X₀) in cm.
    # A longer radiation length is better for a particle detector.
    materials = {
        "A. Titanium": 3.6,
        "B. Aluminium": 8.9,
        "C. 316 Stainless Steel": 1.76,
        "D. Copper": 1.43,
        "E. Nickle": 1.42
    }

    # Find the material with the maximum radiation length
    # The max() function with a key is used on the dictionary's items.
    # The key is a lambda function that specifies we should compare items by their value (the X₀).
    best_choice, max_x0 = max(materials.items(), key=lambda item: item[1])

    # --- Output ---
    print("The unique critical parameter for a particle detector material is to minimize its interaction with particles.")
    print("This is achieved by maximizing the material's radiation length (X₀).\n")
    print("Comparing the radiation lengths of the choices:")

    # Print each material and its radiation length
    for material, x0 in materials.items():
        print(f"Material: {material}, Radiation Length (X₀) = {x0} cm")

    print(f"\nConclusion: {best_choice.split('. ')[1]} is the optimum choice because it has the highest radiation length.")


if __name__ == "__main__":
    solve_detector_material_choice()