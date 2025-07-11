import sys

def solve():
    """
    This function identifies the optimal material for a particle detector's
    cooling system by finding the material with the maximum radiation length.
    """
    # Radiation length (X₀) in cm for common materials.
    # A longer radiation length means less interaction with particles, which is desirable.
    materials = {
        'A. Titanium': 3.56,
        'B. Aluminium': 8.89,
        'C. 316 Stainless Steel': 1.76, # Using value for Iron as a close approximation
        'D. Copper': 1.43,
        'E. Nickle': 1.42
    }

    # Find the material with the highest radiation length
    # The `max` function with a `key` argument is used to find the dictionary key
    # corresponding to the maximum value.
    optimal_material_name = max(materials, key=materials.get)
    optimal_material_value = materials[optimal_material_name]

    print("The unique requirement for a particle detector cooling system is to minimize interaction with the particles being measured.")
    print("This is achieved by selecting the material with the longest radiation length (X₀).\n")
    print("Comparing the radiation lengths:")

    # Print each material and its radiation length
    for name, length in materials.items():
        print(f"- {name}: X₀ = {length} cm")

    print(f"\nThe optimal choice is '{optimal_material_name}' because it has the longest radiation length.")
    
    # Extract the letter from the optimal material name for the final answer
    answer_letter = optimal_material_name.split('.')[0]
    sys.stdout.write(f"\n<<<B>>>\n")


solve()
