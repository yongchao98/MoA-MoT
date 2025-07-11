def solve_electrodynamics_problem():
    """
    This function presents the solution for the electric potential and electric field
    outside a conductive sphere in a uniform electric field in a steady state.
    """
    
    # Define symbolic variables for clarity
    E0 = "E_0"
    R = "R"
    sigma1 = "σ_1"
    sigma2 = "σ_2"
    r = "r"
    theta = "θ"
    r_hat = "r_hat"
    theta_hat = "θ_hat"

    # --- Electric Potential (Phi) for r > R ---
    print("Electric Potential outside the sphere (r > R):")
    
    # Define the coefficient for the dipole term in the potential
    # This corresponds to the number multiplying the R^3/r^2 term
    coeff_phi = f"({sigma1} - {sigma2}) / ({sigma1} + 2*{sigma2})"
    
    # Full expression for the potential
    phi_outside = f"-{E0} * ( (1)*r - {coeff_phi} * {R}^3 / {r}^2 ) * cos({theta})"
    
    print(f"Φ(r, θ) = {phi_outside}\n")
    print("The equation consists of two parts:")
    print(f"1. Potential of the uniform field: - (1) * {E0}*r*cos({theta})")
    print(f"2. Potential of the induced dipole: + ({coeff_phi}) * {E0}*{R}^3/{r}^2 * cos({theta})")
    
    print("-" * 60)

    # --- Electric Field (E) for r > R ---
    print("Electric Field outside the sphere (r > R):")
    print("E(r, θ) = E_r * r_hat + E_θ * θ_hat\n")

    # Radial Component (E_r)
    # Define the coefficient for the dipole term in the radial field component
    coeff_Er = f"2*({sigma1} - {sigma2}) / ({sigma1} + 2*{sigma2})"
    E_r_outside = f"{E0} * ( (1) + {coeff_Er} * {R}^3 / {r}^3 ) * cos({theta})"
    print(f"Radial component E_r = {E_r_outside}")

    # Tangential Component (E_theta)
    # Define the coefficient for the dipole term in the tangential field component
    coeff_E_theta = f"({sigma1} - {sigma2}) / ({sigma1} + 2*{sigma2})"
    E_theta_outside = f"-{E0} * ( (1) - {coeff_E_theta} * {R}^3 / {r}^3 ) * sin({theta})"
    print(f"Tangential component E_θ = {E_theta_outside}\n")
    
    print("The numbers in the field equations are:")
    print(f"For E_r: 1 (from uniform field) and the coefficient 2 in the dipole term's numerator.")
    print(f"For E_θ: 1 (from uniform field) and the coefficient 1 in the dipole term's numerator.")

solve_electrodynamics_problem()