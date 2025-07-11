def solve_magnetostatics_problem():
    """
    This function formats and prints the solution to the magnetostatics problem.
    """
    
    # Define the components of the solution for the magnetic field H(r, theta)
    
    # Inside the sphere (0 < r < R)
    # H_in = ( (2 * mu_0) / mu ) * ( K_0 / (1 + (2 * mu_0 / mu)) ) * z_hat
    
    numerator_in_coeff = "(2 * mu_0)"
    denominator_in_mu_part = "mu"
    numerator_in_K_part = "K_0"
    denominator_in_denom_part = "(1 + (2 * mu_0 / mu))"
    vector_part_in = "z_hat"
    
    inside_field_str = f"H_in = ({numerator_in_coeff} / {denominator_in_mu_part}) * ({numerator_in_K_part} / {denominator_in_denom_part}) * {vector_part_in}"
    
    # As requested by the prompt, here is a printout of each numerical coefficient:
    # From inside field expression: 2, 1, 2
    inside_coeffs = ["2", "1", "2"]
    
    # Outside the sphere (R < r < inf)
    # H_out = ( K_0 / (1 + (2 * mu_0 / mu)) ) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    
    numerator_out_coeff = "K_0"
    denominator_out_coeff = "(1 + (2 * mu_0 / mu))"
    r_dependence = "(R^3 / r^3)"
    angular_dependence = "(2 * cos(theta) * r_hat + sin(theta) * theta_hat)"
    
    outside_field_str = f"H_out = ({numerator_out_coeff} / {denominator_out_coeff}) * {r_dependence} * {angular_dependence}"

    # As requested by the prompt, here is a printout of each numerical coefficient:
    # From outside field expression: 1, 2, 3, 3, 2, 1
    outside_coeffs = ["1", "2", "3", "3", "2", "1"]

    # Construct the final answer string, resembling option E
    final_answer = f"""
H(r, theta) =
  CASE 1: 0 < r < R
    H = (({inside_coeffs[0]} * mu_0) / mu) * (K_0 / ({inside_coeffs[1]} + ({inside_coeffs[2]} * mu_0 / mu))) * z_hat
  
  CASE 2: R < r < inf
    H = (K_0 / ({outside_coeffs[0]} + ({outside_coeffs[1]} * mu_0 / mu))) * (R^{outside_coeffs[2]} / r^{outside_coeffs[3]}) * ({outside_coeffs[4]} * cos(theta) * r_hat + {outside_coeffs[5]} * sin(theta) * theta_hat)
"""

    print("The derived magnetic field is:")
    print(final_answer)

solve_magnetostatics_problem()