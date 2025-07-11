def print_sliding_angle_expression():
    """
    This function prints the derived expression for the tangent of the angle theta
    at which the described rod begins to slide.
    """
    
    print("An expression for the angle θ at which the rod begins to slide can be found by analyzing the forces at the point of impending motion.")
    print("The final relationship is expressed in terms of tan(θ).")
    print("-" * 30)
    
    # Define the numbers in the equation for clarity and to meet the output requirement.
    num_12 = 12
    den_12 = 12
    den_1 = 1
    den_2 = 2
    
    # Construct the string representations of the numerator and denominator
    numerator_str = f"μ * (L^2 / {num_12} + l^2)"
    denominator_str = f"(L^2 / {den_12}) + l^2 * ({den_1} + {den_2}*μ)"
    
    print("The final equation for tan(θ) is:")
    print(f"tan(θ) = ( {numerator_str} ) / ( {denominator_str} )")
    
    print("\nWhere:")
    print("  θ = the angle between the rod and the table")
    print("  μ = the coefficient of static friction")
    print("  L = the total length of the rod")
    print("  l = the distance from the rod's center of mass to the pivot point")

    print("\nTherefore, the expression for the angle θ is:")
    print(f"θ = arctan[ ( {numerator_str} ) / ( {denominator_str} ) ]")

print_sliding_angle_expression()