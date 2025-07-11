def solve_and_print():
    """
    This function prints the derived magnetic field expressions for the given problem.
    It doesn't perform calculations, but presents the final analytical result.
    """
    
    # Final expressions for the magnetic field H in the two regions
    
    # H field inside the sphere (0 < r < R)
    H_in_coefficient = "{2 * mu_0 / mu} * K_0 / (1 + {2 * mu_0 / mu})"
    H_in_vector = "z_hat"
    
    # H field outside the sphere (R < r < infinity)
    H_out_coefficient = "K_0 / (1 + {2 * mu_0 / mu})"
    H_out_radial_dependence = "R^3 / r^3"
    H_out_angular_dependence = "(2 * cos(theta) * r_hat + sin(theta) * theta_hat)"
    
    print("The magnetic field H(r, theta) is given by a piecewise function:")
    print("H(r, theta) =")
    print("{")
    # Print the expression for the inside region
    print(f"  ({H_in_coefficient}) * {H_in_vector},  for 0 < r < R")
    # Print the expression for the outside region
    print(f"  ({H_out_coefficient}) * {H_out_radial_dependence} * {H_out_angular_dependence},  for R < r < infinity")
    print("}")
    print("\nThis result matches Answer Choice E.")

solve_and_print()