def solve_magnetic_shell_problem(H0, R1, R2):
    """
    Calculates the required permeability and interior magnetic field for a
    cylindrical shell that does not distort an external magnetic field.

    Args:
        H0 (float): Magnitude of the externally applied uniform magnetic field (in A/m).
        R1 (float): Inner radius of the cylindrical shell (in m).
        R2 (float): Outer radius of the cylindrical shell (in m).
    """
    if not R2 > R1 > 0:
        print("Error: Radii must be positive and R2 must be greater than R1.")
        return

    # --- Theoretical Result 1: Permeability ---
    # As derived from magnetostatics boundary conditions under the constraint
    # of a non-distorted external field, the relative permeability
    # mu_r = mu/mu_0 must satisfy (mu_r)^2 = 1.
    # Excluding the trivial case mu_r = 1, the required relative permeability is -1.
    mu_r = -1.0

    # --- Theoretical Result 2: Interior Magnetic Field ---
    # For the permeability found above, the magnetic field in the interior
    # region (rho < R1) is uniform and given by H_int = H0 * (R2 / R1)^2.
    H_int_mag = H0 * (R2 / R1)**2

    # --- Output the results ---
    print("This problem asks for the properties of a cylindrical magnetic shell "
          "that prevents distortion of an external uniform magnetic field.")
    print("The analytical solution provides the following results:\n")

    print("1. Required Permeability of the Shell Material")
    print("-----------------------------------------------")
    print("The required permeability 'mu' is related to the vacuum permeability 'mu_0' by:")
    print(f"    mu = mu_r * mu_0")
    print(f"where the required relative permeability, mu_r, is:")
    print(f"    mu_r = {mu_r}")
    print("\nThis non-trivial solution corresponds to a material with negative permeability, "
          "a concept explored in metamaterials for magnetic cloaking.\n")

    print("2. Magnetic Field in the Interior Region (rho < R1)")
    print("-----------------------------------------------------")
    print("The resulting magnetic field inside the shell is uniform, parallel to the "
          "applied field (x-direction), and its magnitude H_int is given by the equation:")
    print(f"    H_int = H0 * (R2 / R1)^2")
    print("\nPlugging in the given values:")
    print(f"    H_int = {H0:.2f} * ({R2:.2f} / {R1:.2f})^2")
    print(f"    H_int = {H0:.2f} * ({R2/R1:.2f})^2")
    print(f"    H_int = {H_int_mag:.2f} A/m")
    print("\nThe final expression for the magnetic field vector in the interior region is:")
    print(f"    H_int = {H_int_mag:.2f} * x_hat (A/m)")


# --- Main execution with example values ---
if __name__ == '__main__':
    # Define example values for the physical parameters
    H0_applied = 10.0   # Magnitude of applied field in A/m
    R1_inner = 0.5      # Inner radius in meters
    R2_outer = 1.0      # Outer radius in meters

    solve_magnetic_shell_problem(H0_applied, R1_inner, R2_outer)