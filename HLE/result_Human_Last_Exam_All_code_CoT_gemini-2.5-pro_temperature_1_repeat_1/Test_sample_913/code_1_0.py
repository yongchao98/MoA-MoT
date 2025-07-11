def print_electric_field_solution():
    """
    This function prints the derived formulas for the electric field
    inside and outside the polarized sensor.
    """
    print("This script presents the analytical solution for the electric field.")
    print("The solution is derived by solving Laplace's equation for the electric potential V")
    print("with boundary conditions determined by the system's geometry and materials.\n")
    
    # --- Region 1: r < R_p ---
    print("--- Electric Field for r < R_p (inside the sensor) ---")
    
    # The electric field vector E is E_r * r_hat + E_theta * theta_hat
    # E = - (A1) * (cos(theta) * r_hat - sin(theta) * theta_hat)
    # where A1 = (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3)
    # The term (cos(theta) * r_hat - sin(theta) * theta_hat) is the unit vector z_hat.
    
    print("The electric field vector is given by:")
    print("E_in = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta) r_hat - sin(theta) theta_hat)\n")
    print("The numerical coefficients and terms in the equation are:")
    print("Coefficient for the entire expression: -1")
    print("Divisor: 3")
    print("Factor: (1 - (R_p**3 / R**3))")
    
    # --- Region 2: R_p < r < R ---
    print("\n--- Electric Field for R_p < r < R (in free space) ---")

    # E_out = - (C1 * r + D1 * r**-2) * (cos(theta)*r_hat - sin(theta)*theta_hat)
    #       - (-2*D1*r**-3*cos(theta)*r_hat + D1*r**-3*sin(theta)*theta_hat)
    # After simplification, the field can be expressed as the sum of a uniform field and a dipole field.
    
    uniform_field_part = "(P_0 / (3 * epsilon_0)) * (R_p/R)**3 * (cos(theta) r_hat - sin(theta) theta_hat)"
    dipole_field_part = "(P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) r_hat + sin(theta) theta_hat)"

    print("The electric field is the superposition of two fields:")
    print("E_out = (Uniform Field from Shell) + (Dipole Field from Sensor)")
    print("\nThe full expression is:")
    print(f"E_out = {uniform_field_part} + {dipole_field_part}\n")
    
    print("The numerical coefficients and terms in the equation are:")
    print("First term (uniform field part):")
    print("  Divisor: 3")
    print("  Factor: (R_p**3 / R**3)")
    print("Second term (dipole part):")
    print("  Divisor: 3")
    print("  Factor: (R_p**3 / r**3)")
    print("  Radial component multiplier: 2")
    print("  Angular component multiplier: 1")

print_electric_field_solution()