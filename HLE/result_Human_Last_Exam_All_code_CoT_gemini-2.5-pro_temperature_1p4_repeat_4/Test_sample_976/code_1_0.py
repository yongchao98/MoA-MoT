def print_solution():
    """
    This function prints the final expressions for the electric potential and
    electric field in the region outside the sphere (r > R).
    """

    # The expressions are derived from first principles and match option B.
    
    potential_outside = "-E_0 * (r - ( (sigma_1 - sigma_2) * R**3 ) / ( (sigma_1 + 2*sigma_2) * r**2 )) * cos(theta)"
    
    # Electric field is a vector with r and theta components
    E_r_component = "E_0 * (1 + ( 2*(sigma_1 - sigma_2) * R**3 ) / ( (sigma_1 + 2*sigma_2) * r**3 )) * cos(theta) * r_hat"
    E_theta_component = "- E_0 * (1 - ( (sigma_1 - sigma_2) * R**3 ) / ( (sigma_1 + 2*sigma_2) * r**3 )) * sin(theta) * theta_hat"

    print("The electric potential outside the sphere (r > R) is:")
    print(f"Phi(r, theta) = {potential_outside}\n")
    
    print("The electric field outside the sphere (r > R) is:")
    print(f"E(r, theta) = {E_r_component} {E_theta_component}")

print_solution()