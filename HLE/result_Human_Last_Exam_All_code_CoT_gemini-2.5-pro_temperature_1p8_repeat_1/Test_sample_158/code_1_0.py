import sympy

def solve_magnetic_shielding():
    """
    This function derives and prints the symbolic solutions for the permeability
    of a cylindrical shell and the resulting interior magnetic field required
    to prevent distortion of an external uniform magnetic field.

    The derivation is based on solving Laplace's equation for the magnetic
    scalar potential and applying boundary conditions, under the key insight that
    the interior region (rho < R1) must be a perfect diamagnet (mu_int = 0)
    to achieve a non-trivial solution.
    """

    # Define symbolic variables for clarity
    R1, R2, H0, mu0 = sympy.symbols("R1 R2 H0 mu0")

    print("This script provides the symbolic solution to the magnetic shielding problem.\n")

    # --- Part 1: Permeability of the Shell ---
    print("1. Required Permeability of the Shell (mu)")
    print("-" * 40)
    print("For the external field to remain undistorted, and excluding the trivial")
    print("case (mu = mu0), the shell must have a specific permeability `mu` that")
    print("depends on the inner and outer radii, R1 and R2.")
    print("The derived value, assuming a perfect diamagnetic core (mu_int = 0), is:\n")

    # Construct and print the permeability equation
    mu_numerator = R2**2 + R1**2
    mu_denominator = R2**2 - R1**2
    mu_expr = mu0 * (mu_numerator / mu_denominator)

    print(f"The final equation for the permeability mu is:")
    print(f"mu = {mu0} * ({mu_numerator}) / ({mu_denominator})\n")

    # --- Part 2: Interior Magnetic Field ---
    print("2. Magnetic Field in the Interior (H_int)")
    print("-" * 40)
    print("Under this condition, the magnetic field `H_int` inside the shell (rho < R1)")
    print("is uniform and aligned with the external field (in the x-direction).")
    print("Its magnitude is given by:\n")

    # Construct and print the interior field equation
    h_int_numerator = 2 * H0 * R2**2
    h_int_denominator = R1**2 + R2**2
    
    print(f"The final equation for the interior field H_int is:")
    print(f"H_int = [({h_int_numerator}) / ({h_int_denominator})] * x_hat\n")

    print("Let's break down the numerical factors in the H_int equation:")
    print(f"  - The numerator of the coefficient is: {h_int_numerator}")
    print(f"  - The denominator of the coefficient is: {h_int_denominator}")
    print("  - The constant '2' appears explicitly in the numerator of the coefficient.")

if __name__ == '__main__':
    solve_magnetic_shielding()