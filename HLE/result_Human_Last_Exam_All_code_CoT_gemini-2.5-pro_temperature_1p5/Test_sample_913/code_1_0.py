def print_electric_fields():
    """
    This function prints the derived expressions for the electric field
    in the two specified regions.
    """

    # --- Region 1: Inside the sensor (r < R_p) ---
    print("Electric field for r < R_p (inside the sensor):")
    # The term (cos(theta) r_hat - sin(theta) theta_hat) is the unit vector z_hat.
    # The field is uniform and points in the -z direction.
    E_in_str = (
        "E_in = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * "
        "(cos(theta) r_hat - sin(theta) theta_hat)"
    )
    print(E_in_str)
    print("\n" + "="*50 + "\n")


    # --- Region 2: Between sensor and shell (R_p < r < R) ---
    print("Electric field for R_p < r < R (in free space):")
    # The field is a superposition of a uniform field and a dipole field.
    E_out_uniform_part = (
        "(P_0 / (3 * epsilon_0)) * (R_p/R)**3 * "
        "(cos(theta) r_hat - sin(theta) theta_hat)"
    )
    E_out_dipole_part = (
        "(P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * "
        "(2*cos(theta) r_hat + sin(theta) theta_hat)"
    )
    E_out_str = f"E_out = {E_out_uniform_part} + {E_out_dipole_part}"
    print(E_out_str)


if __name__ == "__main__":
    print_electric_fields()
    # Comparing our rigorously derived fields with the options provided,
    # Choice E correctly describes the field inside the sensor (r < R_p).
    # The other choices have a sign error for this region, making them
    # physically incorrect. Although Choice E contains a sign error in the
    # expression for the field outside the sensor, it is the best fit.
    # Final Answer is E.
