def solve_magnetic_shell_problem():
    """
    This function presents the solution to the specified magnetostatics problem.
    It calculates and prints the required permeability of a cylindrical shell
    and the resulting magnetic field in its interior, under the condition that
    the external magnetic field is not distorted.
    """

    # --- Introduction to the Solution ---
    print("This script provides the solution to the magnetic cylindrical shell problem.")
    print("The goal is to find the shell's permeability (mu) and the interior magnetic field (H_int).")
    print("The key condition is that the external magnetic field H_ext remains H0*x_hat.\n")

    # --- Part 1: Required Permeability ---
    print("--------------------------------------------------")
    print("1. Required Permeability of the Shell Material")
    print("--------------------------------------------------")
    print("Analysis of the boundary conditions reveals that a non-trivial solution")
    print("(where mu is not mu_0 and the shell has finite thickness) requires")
    print("a specific value for the permeability of the material in the shell.")
    print("\nThe required permeability, mu, is found to be:")
    print("  mu = -1 * mu_0")
    print("\nWhere mu_0 is the permeability of free space.")
    print("The numerical constant in this equation is -1.")
    print("This corresponds to a relative permeability mu_r = -1, which is characteristic")
    print("of certain exotic materials or metamaterials.\n")

    # --- Part 2: Magnetic Field in the Interior Region ---
    print("--------------------------------------------------")
    print("2. Magnetic Field in the Interior Region (rho < R1)")
    print("--------------------------------------------------")
    print("Using the determined permeability (mu = -mu_0), the magnetic field")
    print("in the interior region, H_int, is found to be uniform and is given by:")
    print("  H_int = H0 * (R2 / R1)**2 * x_hat")
    print("\nWhere:")
    print("  H0 is the magnitude of the externally applied uniform magnetic field.")
    print("  R1 is the interior radius of the shell.")
    print("  R2 is the exterior radius of the shell.")
    print("  x_hat is the unit vector in the x-direction.")
    print("\nThe numerical constants in this equation are:")
    print("  - The exponent '2' applied to the ratio of the radii.")
    # The prompt asks to output each number in the equation. While R1 and R2 are symbols,
    # the indices 1 and 2 are numbers that identify them.
    print("  - The index '2' identifying the exterior radius R2.")
    print("  - The index '1' identifying the interior radius R1.")

if __name__ == '__main__':
    solve_magnetic_shell_problem()