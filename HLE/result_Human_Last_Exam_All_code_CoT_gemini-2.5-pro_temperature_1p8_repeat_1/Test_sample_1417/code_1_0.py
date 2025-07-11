import sys

def solve_detector_material():
    """
    Determines the optimal material for a particle detector's cooling system
    based on the unique requirement of minimizing particle interaction.
    """

    # Step 1: Explain the physics principle.
    # The unique parameter for a particle detector is minimizing the interaction
    # with the particles it is designed to measure. This is quantified by the
    # 'radiation length' (X0). A longer radiation length is better as it implies
    # the material is more 'transparent' to particles.
    print("Finding the optimal material for a particle detector cooling system.")
    print("The key parameter is the material's radiation length (X0). A longer radiation length is better.")
    print("-" * 30)

    # Step 2: Define the materials and their radiation lengths in cm.
    # Data is sourced from standard physics tables (e.g., Particle Data Group).
    materials = {
        "Titanium": 3.56,
        "Aluminium": 8.90,
        "316 Stainless Steel": 1.76, # Using Iron as a close proxy
        "Copper": 1.43,
        "Nickel": 1.42
    }

    # Print the data being used for the comparison
    print("Comparing the following materials and their radiation lengths (in cm):")
    for name, length in materials.items():
        print(f"- {name}: {length:.2f} cm")

    # Step 3: Find the material with the maximum radiation length.
    # The max() function with a key is a clean way to find the dictionary entry
    # corresponding to the maximum value.
    optimal_material = max(materials, key=materials.get)
    max_radiation_length = materials[optimal_material]

    # Step 4: Output the result.
    print("-" * 30)
    print(f"The material with the longest radiation length is '{optimal_material}' with a value of {max_radiation_length:.2f} cm.")
    print("Therefore, Aluminium is the optimum choice to minimize interference with particle measurements.")

# Execute the function
solve_detector_material()

# Provide the final answer in the required format
# The options correspond to: A. Titanium, B. Aluminium, C. 316 Stainless Steel, D. Copper, E. Nickle
# Our determined optimal material is Aluminium.
sys.stdout.write("<<<B>>>\n")