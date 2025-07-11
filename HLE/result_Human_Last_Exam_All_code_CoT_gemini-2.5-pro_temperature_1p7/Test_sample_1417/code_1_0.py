def find_optimal_material():
    """
    Identifies the optimal material for a particle detector's cooling tubes
    by maximizing the radiation length (X0) to ensure transparency to particles.
    """
    # Data sourced from the Particle Data Group (PDG)
    # The unique parameter to optimize is the radiation length (X0).
    # A longer radiation length means the material is more 'transparent' to particles,
    # causing less interference with detector measurements.
    materials_data = {
        'A. Titanium': 3.56,         # Radiation Length in cm
        'B. Aluminium': 8.9,          # Radiation Length in cm
        'C. 316 Stainless Steel': 1.76, # Approx. value for Iron, in cm
        'D. Copper': 1.43,          # Radiation Length in cm
        'E. Nickel': 1.48           # Radiation Length in cm
    }

    print("The unique parameter for a particle detector is minimizing particle interaction.")
    print("This is achieved by maximizing the radiation length (X0).\n")
    print("Comparing the radiation lengths for the given materials:")
    
    # Initialize variables to find the maximum
    optimal_material_name = None
    max_x0 = -1

    # Loop through the data to find the max and print each value
    for material, x0 in materials_data.items():
        print(f"Material: {material}, Radiation Length (X0) = {x0} cm")
        if x0 > max_x0:
            max_x0 = x0
            optimal_material_name = material
    
    # Extract the choice letter from the name
    optimal_choice_letter = optimal_material_name.split('.')[0]

    print(f"\nThe material with the maximum radiation length is {optimal_material_name} with X0 = {max_x0} cm.")
    print(f"Therefore, Aluminium is the optimum choice.")
    
    # Final answer in the required format
    print(f"\n<<<{optimal_choice_letter}>>>")

find_optimal_material()