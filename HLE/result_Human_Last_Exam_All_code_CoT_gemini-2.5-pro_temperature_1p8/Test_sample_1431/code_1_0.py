def solve_rope_problem():
    """
    This function prints the derived formula for the force F required to lift the rope.
    
    The formula is derived from the work-energy theorem:
    Work done by F = Change in Kinetic Energy + Change in Potential Energy
    F * l = (1/2 * m * v^2) + (m * g * l / 2)
    
    Solving for F gives:
    F = (m * v^2) / (2 * l) + (m * g) / 2
    """
    
    # Define the variables symbolically for the printout
    m = "m"  # mass of the rope
    l = "l"  # length of the rope
    v = "v"  # final speed of the rope
    g = "g"  # acceleration due to gravity
    
    # The numbers in the final equation are 2.
    denominator = 2
    
    # Print the final derived equation for the force F
    print("The mystical challenge is solved! The formula for the required force F is:")
    print(f"F = ({m} * {v}^2) / ({denominator} * {l}) + ({m} * {g}) / {denominator}")
    print("\nWhere:")
    print("F = The constant force applied by the monk")
    print("m = Mass of the rope")
    print("l = Length of the rope")
    print("v = The speed of the rope as it just leaves the ground")
    print("g = The acceleration due to gravity (approximately 9.8 m/s^2)")

# Execute the function to display the answer
solve_rope_problem()