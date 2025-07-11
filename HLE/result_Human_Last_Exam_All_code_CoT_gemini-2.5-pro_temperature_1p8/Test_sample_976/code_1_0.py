def solve_electrodynamics_problem():
    """
    This function prints the derived expressions for the electric potential and 
    electric field outside the sphere for the given problem.
    """

    # --- Symbolic Representation of Variables ---
    E0 = "E_0"
    sigma1 = "sigma_1"
    sigma2 = "sigma_2"
    R = "R"
    r = "r"
    theta = "theta"
    r_hat = "r_hat"  # unit vector in radial direction
    theta_hat = "theta_hat" # unit vector in polar angle direction
    
    # --- Introduction ---
    print("This script provides the solution for the electric potential and field in the region outside the sphere (r > R).")
    print("The derivation is based on solving Laplace's equation with steady-state boundary conditions at the sphere's surface.")
    print("-" * 70)

    # --- Electric Potential (Phi) for r > R ---
    # The potential is the sum of the external field potential and a dipole term.
    # Phi_out = -E0 * r * cos(theta) + (dipole_term) * cos(theta)
    phi_dipole_coeff = f"({sigma1} - {sigma2}) * {R}^3 / (({sigma1} + 2*{sigma2})*{r}^2)"
    phi_out_expression = f"-{E0} * (r - ({sigma1} - {sigma2}) * {R}^3 / (({sigma1} + 2*{sigma2}) * {r}^2)) * cos({theta})"

    print("\nElectric Potential Phi(r, theta) for r > R:")
    # Print each term of the equation explicitly
    print(f"Phi(r > R) = -{E0} * (r - (({sigma1} - {sigma2}) * {R}^3) / (({sigma1} + 2 * {sigma2}) * {r}^2)) * cos({theta})")


    # --- Electric Field (E) for r > R ---
    # E = -grad(Phi). It has a radial (E_r) and a theta (E_theta) component.
    
    # E_r component
    Er_factor_numerator = f"2*({sigma1} - {sigma2})*{R}^3"
    Er_factor_denominator = f"({sigma1} + 2*{sigma2})*{r}^3"
    Er_expression = f"{E0} * (1 + {Er_factor_numerator} / ({Er_factor_denominator})) * cos({theta}) {r_hat}"

    # E_theta component
    Etheta_factor_numerator = f"({sigma1} - {sigma2})*{R}^3"
    Etheta_factor_denominator = f"({sigma1} + 2*{sigma2})*{r}^3"
    Etheta_expression = f"-{E0} * (1 - {Etheta_factor_numerator} / ({Etheta_factor_denominator})) * sin({theta}) {theta_hat}"
    
    print("\n\nElectric Field E(r, theta) for r > R:")
    # Print the radial component of the field
    print(f"E_r = {E0} * (1 + (2 * ({sigma1} - {sigma2}) * {R}^3) / (({sigma1} + 2 * {sigma2}) * {r}^3)) * cos({theta}) * {r_hat}")
    # Print the theta component of the field
    print(f"E_theta = -{E0} * (1 - (({sigma1} - {sigma2}) * {R}^3) / (({sigma1} + 2 * {sigma2}) * {r}^3)) * sin({theta}) * {theta_hat}")
    # Print the full vector expression
    print(f"\nE_vector(r > R) = [ {E0} * (1 + (2 * ({sigma1} - {sigma2}) * {R}^3) / (({sigma1} + 2 * {sigma2}) * {r}^3)) * cos({theta}) ] {r_hat} \n                 - [ {E0} * (1 - (({sigma1} - {sigma2}) * {R}^3) / (({sigma1} + 2 * {sigma2}) * {r}^3)) * sin({theta}) ] {theta_hat}")


solve_electrodynamics_problem()