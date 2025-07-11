def solve_electrodynamics_problem():
    """
    This function prints the expressions for the electric potential and electric field
    in the region outside the sphere (r > R).
    """

    # Define the components of the equations for clarity.
    # The common factor in the coefficients
    coeff_factor = "(sigma_1 - sigma_2) / (sigma_1 + 2*sigma_2)"
    
    # Expression for the electric potential Phi(r, theta) for r > R
    potential_out = f"-E_0 * (r - ({coeff_factor}) * R^3/r^2) * cos(theta)"
    
    # Expression for the radial component of the Electric Field E_r(r, theta) for r > R
    e_field_r_coeff = f"1 + 2 * ({coeff_factor}) * R^3/r^3"
    e_field_r = f"E_0 * ({e_field_r_coeff}) * cos(theta) * r_hat"

    # Expression for the angular component of the Electric Field E_theta(r, theta) for r > R
    e_field_theta_coeff = f"1 - ({coeff_factor}) * R^3/r^3"
    e_field_theta = f"-E_0 * ({e_field_theta_coeff}) * sin(theta) * theta_hat"

    print("The electric potential and electric field in the region outside the sphere (r > R) are:")
    print("-" * 80)
    
    # Print the final equations, showing each part as requested
    print("Potential Phi(r, theta) for r > R:")
    print(f"Phi(r, theta) = {potential_out}")
    print("\n")
    
    print("Electric Field E(r, theta) for r > R:")
    print("E(r, theta) = E_r + E_theta")
    print(f"E_r = {e_field_r}")
    print(f"E_theta = {e_field_theta}")
    print("-" * 80)
    print("\nThese expressions correspond to Option B.")

solve_electrodynamics_problem()