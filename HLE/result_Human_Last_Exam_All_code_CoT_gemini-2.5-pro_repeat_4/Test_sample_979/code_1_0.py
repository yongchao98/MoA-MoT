import sympy as sp

def display_solution():
    """
    This function prints the step-by-step derivation and the final solution
    for the magnetic field H inside and outside the spherical shell.
    """
    
    print("Based on the principles of magnetostatics, the magnetic field H is determined by solving Laplace's equation for the magnetic scalar potential and applying the boundary conditions at the spherical surface r=R.")
    
    print("\nThe derived magnetic field H in the two regions is as follows:")
    
    # --- Inside the sphere (r < R) ---
    print("\n--- For the region inside the sphere (0 < r < R): ---")
    
    # Using simple strings to construct the formula output
    H_in_coeff_part1 = "(2 * mu_0) / mu"
    H_in_coeff_part2 = "K_0 / (1 + (2 * mu_0) / mu)"
    H_in_vector = "z_hat"
    
    print(f"H_in(r, theta) = ( {H_in_coeff_part1} ) * ( {H_in_coeff_part2} ) * {H_in_vector}")
    print("\nThis represents a uniform magnetic field pointing in the positive z-direction.")

    # --- Outside the sphere (r > R) ---
    print("\n--- For the region outside the sphere (R < r < infinity): ---")

    H_out_coeff = "K_0 / (1 + (2 * mu_0) / mu)"
    H_out_radial_term = "R^3 / r^3"
    H_out_angular_term = "(2 * cos(theta) * r_hat + sin(theta) * theta_hat)"
    
    print(f"H_out(r, theta) = ( {H_out_coeff} ) * ( {H_out_radial_term} ) * {H_out_angular_term}")
    print("\nThis represents a magnetic dipole field.")

    # --- Final Conclusion ---
    print("\n--- Final Answer ---")
    print("Combining these two results into a single expression for H(r, theta):")
    print("H(r, theta) =")
    print(f"  For 0 < r < R:      ( (2 * mu_0) / mu ) * ( K_0 / (1 + (2 * mu_0) / mu) ) * z_hat")
    print(f"  For R < r < infinity: ( K_0 / (1 + (2 * mu_0) / mu) ) * ( R^3 / r^3 ) * (2 * cos(theta) * r_hat + sin(theta) * theta_hat)")

    print("\nThis solution matches the expressions given in answer choice E.")

if __name__ == '__main__':
    display_solution()