def calculate_primitives_for_methane():
    """
    Calculates the total number of primitive Gaussians for a methane (CH4)
    molecule using the 6-311G** basis set.
    """
    
    # Define the number of primitives for each type of atom based on 6-311G**
    
    # For Hydrogen (H):
    # s-primitives (from '311' valence): 3 + 1 + 1 = 5
    # p-primitives (from '**'): 3
    primitives_H = 5 + 3
    
    # For a second-row heavy atom like Carbon (C):
    # s-primitives (1s core '6' + 2s valence '311'): 6 + 3 + 1 + 1 = 11
    # p-primitives (2p valence '311'): 3 + 1 + 1 = 5
    # d-primitives (from '*'): 6
    primitives_C = 11 + 5 + 6
    
    # Methane (CH4) has 1 Carbon atom and 4 Hydrogen atoms
    num_C = 1
    num_H = 4
    
    # Calculate total primitives
    total_primitives = (num_C * primitives_C) + (num_H * primitives_H)
    
    # Print the explanation and calculation
    print("Calculation of primitive Gaussians for Methane (CH4) with the 6-311G** basis set:")
    print("-" * 75)
    print(f"Primitives for one Carbon (C) atom: {primitives_C}")
    print(f"Primitives for one Hydrogen (H) atom: {primitives_H}")
    print("-" * 75)
    
    # Build and print the equation for methane
    equation_parts = []
    # Add carbon part
    equation_parts.append(str(primitives_C) + " (C)")
    # Add hydrogen parts
    for _ in range(num_H):
        equation_parts.append(str(primitives_H) + " (H)")
        
    equation_str = " + ".join(equation_parts)
    print(f"Total for CH4 = {equation_str}")
    print(f"Total = {primitives_C} + {num_H * primitives_H}")
    print(f"Final total number of primitive Gaussians = {total_primitives}")

# Run the calculation
calculate_primitives_for_methane()