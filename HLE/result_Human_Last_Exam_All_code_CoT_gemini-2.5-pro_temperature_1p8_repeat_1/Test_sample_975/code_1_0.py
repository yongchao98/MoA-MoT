def display_final_equations():
    """
    This function prints the final derived expressions for the magnetic field H
    in both regions, as requested by the user.
    """
    
    # Header for the solution
    print("Based on the physical derivation, the magnetic field H(r, θ) is as follows:\n")
    
    # Region 1: Inside the shield (0 < r < Rp)
    print("="*60)
    print("In the region 0 < r < R_p (inside the shield):")
    print("="*60)
    
    # The problem asks to output the numbers in the equation.
    # The overall coefficient is M_0 * (2*R_p^3 + 1*R^3) / (3*R^3)
    # The vector part is (-1*cos(theta)*i_r + 1*sin(theta)*i_theta)
    
    print("H = M_0 * (2*R_p^3 + R^3) / (3*R^3) * (-cos(theta) * i_r + sin(theta) * i_theta)")
    print("\nThe numerical coefficients in this expression are 2, 1, and 3.\n")
    
    # Region 2: Between the shield and the conductor (Rp < r < R)
    print("="*60)
    print("In the region R_p < r < R (between the shield and conductor):")
    print("="*60)
    
    # The problem asks to output the numbers in the equation.
    # The radial component coefficient is -2/3.
    # The theta component coefficient is 1/3, with an internal factor of 2.
    
    print("H_r = (-2/3) * M_0 * [ (R_p/R)^3 - (R_p/r)^3 ] * cos(theta)")
    print("H_theta = (1/3) * M_0 * [ 2*(R_p/R)^3 + (R_p/r)^3 ] * sin(theta)")
    print("\nThe numerical coefficients in the components are -2/3, 1, 1/3, and 2.")
    print("\nThese results correspond to Answer Choice B.")

# Execute the function to print the solution
display_final_equations()