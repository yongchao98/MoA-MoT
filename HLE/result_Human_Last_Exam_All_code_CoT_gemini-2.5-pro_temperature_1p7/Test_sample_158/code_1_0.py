import sympy

def solve_magnetic_shell_problem():
    """
    This function presents the analytical solution to the magnetic cylindrical shell problem.

    The problem asks for the permeability (mu) of a cylindrical shell (radii R1, R2)
    such that it produces no distortion to an external uniform magnetic field (H0)
    in the exterior region, and to find the resulting magnetic field (H_int)
    inside the shell.
    """

    # Define symbolic variables to represent the physical quantities for clarity
    H0_mag, R1, R2, mu0 = sympy.symbols('H0 R1 R2 mu_0')
    
    # --- Part 1: Required Permeability ---
    # The derivation by applying magnetostatic boundary conditions leads to the
    # following algebraic constraint for the relative permeability mu_r = mu / mu_0:
    # (mu_r**2 - 1) * (R1**2 - R2**2) = 0
    #
    # Since we have a shell, R1 is not equal to R2, so (R1**2 - R2**2) is non-zero.
    # Therefore, we must have (mu_r**2 - 1) = 0.
    # This equation has two solutions: mu_r = 1 (trivial case) and mu_r = -1.
    # The problem asks for the non-trivial solution.

    permeability_mu_r = -1

    print("--- Required Permeability of the Shell Material ---")
    print("For the external magnetic field to be undistorted, the permeability 'mu' of the shell must be:")
    print(f"μ = {permeability_mu_r} * μ₀")
    print("(where μ₀ is the permeability of free space)")

    print("\n" + "-"*50 + "\n")

    # --- Part 2: Magnetic Field in the Interior Region ---
    # With the determined value of mu_r = -1, the system of boundary
    # condition equations can be solved for the magnetic field inside the cavity.
    # The resulting field H_int is uniform and aligned with the external field H0.
    # Its magnitude is found to be amplified.

    print("--- Magnetic Field in the Interior Region (ρ < R1) ---")
    print("The magnetic field vector 'H_int' inside the shell is given by the equation:")
    print("H_int = H₀ * (R₂² / R₁²) * x̂")
    print("\nBreaking down the equation for the magnitude of the interior field:")
    print("Magnitude(H_int) = H₀ * (R₂² / R₁²)")
    print("\nIn this equation:")
    print("- H₀ is the magnitude of the externally applied uniform magnetic field.")
    print("- R₁ is the interior radius of the shell.")
    print("- R₂ is the exterior radius of the shell.")
    print("- The number '2' is the exponent for the radii.")

if __name__ == '__main__':
    solve_magnetic_shell_problem()
