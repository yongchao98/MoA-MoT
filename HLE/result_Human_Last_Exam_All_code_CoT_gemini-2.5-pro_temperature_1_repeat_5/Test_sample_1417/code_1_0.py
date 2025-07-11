import operator

def solve_detector_material_choice():
    """
    Determines the optimal material for a particle detector's cooling system
    based on the unique parameter of minimizing particle interaction.
    """
    print("Analyzing the optimal material choice for the particle detector cooling system.\n")

    # Step 1 & 2: Explain the reasoning and the unique parameter
    print("Step 1: Identify the unique parameter.")
    print("The primary function of a particle detector is to measure particles accurately.")
    print("Any material in the particle's path, such as the cooling tubes, will interfere with this measurement.")
    print("The unique requirement for this application, compared to a spacecraft, is to minimize this interference.\n")

    print("Step 2: Quantify the parameter.")
    print("This property is quantified by the material's 'radiation length' (X0).")
    print("A longer radiation length means particles can travel through it with less interaction.")
    print("Therefore, the optimum material is the one with the maximum radiation length.\n")

    # Step 3: Define material data
    # Radiation length (X0) values in cm from the Particle Data Group (PDG).
    # For 316 Stainless Steel, the value is approximated based on its main component, Iron (X0 = 1.76 cm).
    materials_data = {
        "A. Titanium": 3.56,
        "B. Aluminium": 8.9,
        "C. 316 Stainless Steel": 1.74,
        "D. Copper": 1.43,
        "E. Nickle": 1.42
    }

    print("Step 3: Compare the radiation lengths (X0) of the candidate materials.")
    print("The material with the highest X0 value is the best choice.")
    print("-" * 40)
    print("Material              | Radiation Length (X0) in cm")
    print("-" * 40)
    for material, x0 in materials_data.items():
        print(f"{material:<22}| {x0}")
    print("-" * 40)
    print("\n")


    # Step 4: Find the optimum material
    # The key in the dictionary already contains the letter designation.
    optimum_material_name, optimum_x0 = max(materials_data.items(), key=operator.itemgetter(1))

    print("Step 4: Conclusion.")
    print(f"Based on the data, {optimum_material_name.split('. ')[1]} has the highest radiation length of {optimum_x0} cm.")
    print("This means it will interfere the least with the particles being measured.")
    print(f"Therefore, considering only this unique parameter, {optimum_material_name.split('. ')[1]} is the optimum choice.")


solve_detector_material_choice()
# The final answer is the letter corresponding to the best material.
# From the analysis, Aluminium (B) has the highest radiation length.
print("<<<B>>>")