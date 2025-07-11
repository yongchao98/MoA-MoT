def print_solution():
    """
    This function prints the final derived equations for the magnetic field H.
    The expressions match choice B from the problem description.
    """
    
    print("The magnetic field H in the two specified regions is determined as follows:")
    print("-" * 70)
    
    # Expression for the field inside the shield
    print("In the region 0 < r < R_p:")
    # We break down the expression to show all numerical components as requested.
    # The term (-cos(theta) i_r + sin(theta) i_theta) is equivalent to -z_hat.
    # The coefficient is M_0 * (2*R_p^3 + R^3) / (3*R^3)
    print("  H = M_0 * ( (2*R_p^3 + R^3) / (3*R^3) ) * ( - cos(theta) * i_r + sin(theta) * i_theta )")
    
    print("\n") # Add a newline for better readability
    
    # Expression for the field between the shield and the conductor
    print("In the region R_p < r < R:")
    # The field is broken into its radial (r) and polar (theta) components.
    # H_r = - (2*M_0 / 3) * [ (R_p/R)^3 - (R_p/r)^3 ] * cos(theta)
    # H_theta = (M_0 / 3) * [ 2*(R_p/R)^3 + (R_p/r)^3 ] * sin(theta)
    print("  H = - (2*M_0 / 3) * [ (R_p/R)^3 - (R_p/r)^3 ] * cos(theta) * i_r \\")
    print("      + (M_0 / 3) * [ 2*(R_p/R)^3 + (R_p/r)^3 ] * sin(theta) * i_theta")
    
    print("-" * 70)
    print("\nThese expressions correspond to Answer Choice B.")

# Execute the function to print the solution
print_solution()