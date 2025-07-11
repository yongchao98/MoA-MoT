def solve_electrostatics_problem():
    """
    This function presents the solution for the electric potential and electric field
    outside a conductive sphere in a uniform electric field in a steady state.
    """

    # The problem asks for the potential and field for r > R.
    # The following strings represent the components of the correct formulas.
    
    # --- Electric Potential Phi(r, theta) for r > R ---
    print("The electric potential Phi(r, theta) for r > R is:")
    
    # The potential outside is composed of the applied uniform field potential and a dipole term.
    # Phi_out = -E0 * (r - (sigma1 - sigma2)*R^3 / ((sigma1 + 2*sigma2)*r^2)) * cos(theta)
    
    term_uniform_potential = "-E_0 * r * cos(theta)"
    
    # The coefficient of the dipole term
    dipole_coefficient_numerator = "(sigma_1 - sigma_2) * R^3"
    dipole_coefficient_denominator = "(sigma_1 + 2*sigma_2) * r^2"
    
    # Full expression for printing
    phi_expression = f"Phi(r > R, theta) = -E_0 * (r - ({dipole_coefficient_numerator}) / ({dipole_coefficient_denominator})) * cos(theta)"
    print(phi_expression)
    print("\n")

    # --- Electric Field E(r, theta) for r > R ---
    print("The electric field E(r, theta) for r > R is:")
    
    # The field is the negative gradient of the potential, E = -grad(Phi).
    # It has a radial (r_hat) and an angular (theta_hat) component.

    # Radial component E_r
    # E_r = E0 * (1 + 2*(sigma1 - sigma2)*R^3 / ((sigma1 + 2*sigma2)*r^3)) * cos(theta)
    E_r_dipole_numerator = "2 * (sigma_1 - sigma_2) * R^3"
    E_r_dipole_denominator = "(sigma_1 + 2*sigma_2) * r^3"
    E_r_expression = f"E_r = E_0 * (1 + ({E_r_dipole_numerator}) / ({E_r_dipole_denominator})) * cos(theta)"

    # Angular component E_theta
    # E_theta = -E0 * (1 - (sigma1 - sigma2)*R^3 / ((sigma1 + 2*sigma2)*r^3)) * sin(theta)
    E_theta_dipole_numerator = "(sigma_1 - sigma_2) * R^3"
    E_theta_dipole_denominator = "(sigma_1 + 2*sigma_2) * r^3"
    E_theta_expression = f"E_theta = -E_0 * (1 - ({E_theta_dipole_numerator}) / ({E_theta_dipole_denominator})) * sin(theta)"

    # Full vector expression
    E_vector_expression = f"E(r > R, theta) = ({E_r_expression}) r_hat + ({E_theta_expression}) theta_hat"
    
    print("The radial component is:")
    print(f"E_r = E_0 * (1 + 2 * (sigma_1 - sigma_2) * R^3 / ((sigma_1 + 2*sigma_2) * r^3)) * cos(theta)")
    print("\nThe angular component is:")
    print(f"E_theta = -E_0 * (1 - (sigma_1 - sigma_2) * R^3 / ((sigma_1 + 2*sigma_2) * r^3)) * sin(theta)")

solve_electrostatics_problem()