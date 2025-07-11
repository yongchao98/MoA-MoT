import sys

def solve_detector_material():
    """
    Identifies the optimum material for a particle detector's cooling system
    by maximizing the radiation length (X0), which minimizes particle interaction.
    """
    # The unique parameter for a particle detector is minimizing interaction
    # with the particles being measured. This is quantified by the 'radiation length' (X0).
    # A longer radiation length is better, as it means the material is more 'transparent'
    # to the particles. The goal is to maximize this parameter.

    # Data Source: Particle Data Group (PDG). Values are radiation length (X0) in cm.
    # Generally, materials with lower atomic number (Z) and lower density have longer radiation lengths.
    material_properties = {
        "A. Titanium": 3.56,
        "B. Aluminium": 8.9,
        "C. 316 Stainless Steel": 1.74, # Dominated by Iron (Fe)
        "D. Copper": 1.43,
        "E. Nickle": 1.42
    }

    # Initialize variables to find the best material
    optimum_material_name = ""
    max_radiation_length = -1.0

    print("Step 1: Define the unique optimization parameter.")
    print("For a particle detector, the most important unique parameter for internal components is maximizing 'particle transparency'.")
    print("This is achieved by choosing the material with the longest radiation length (X0).\n")

    print("Step 2: Compare the radiation lengths of the candidate materials.")
    # Iterate through the dictionary to find the maximum value
    for material, length in material_properties.items():
        print(f"Material: {material}, Radiation Length (X0) = {length} cm")
        if length > max_radiation_length:
            max_radiation_length = length
            optimum_material_name = material

    print("\nStep 3: Identify the optimum choice.")
    print("The material with the maximum radiation length will interfere the least with particle measurements.")
    print(f"The optimum material is {optimum_material_name}, with the longest radiation length of {max_radiation_length} cm.")

solve_detector_material()
<<<B>>>