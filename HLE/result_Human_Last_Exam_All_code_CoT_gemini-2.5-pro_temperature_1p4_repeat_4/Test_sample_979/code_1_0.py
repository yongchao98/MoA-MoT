def display_magnetic_field_solution():
    """
    This function prints the final derived expression for the magnetic field H(r, theta)
    for both the interior and exterior of the spherical shell.
    """
    
    # The expressions are constructed from the derivation.
    # Note that z_hat is used for the uniform field inside for conciseness.
    H_inside = "    (2 * mu_0 / mu) * (K_0 / (1 + (2 * mu_0 / mu))) * z_hat"
    H_outside = "    (K_0 / (1 + (2 * mu_0 / mu))) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"

    print("The solution for the magnetic field H(r, theta) is:")
    print("H(r, theta) =")
    print("{")
    print(f"  (line 1) {H_inside}      , for 0 < r < R")
    print(f"  (line 2) {H_outside}  , for R < r < infinity")
    print("}")
    print("\nThis result matches answer choice E.")

if __name__ == '__main__':
    display_magnetic_field_solution()