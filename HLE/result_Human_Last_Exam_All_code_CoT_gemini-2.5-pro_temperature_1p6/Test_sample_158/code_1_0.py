def solve_magnetostatics_problem():
    """
    This function presents the solution to the cylindrical shell magnetostatics problem.
    It calculates and prints the required permeability (mu) of the shell and the
    resulting magnetic field (H_int) in the interior region.
    The results are derived analytically and presented as symbolic formulas.
    """

    # --- Symbolic placeholders for the given parameters ---
    mu_0 = "μ₀"  # Permeability of free space
    H_0 = "H₀"    # Magnitude of the externally applied uniform magnetic field
    R_1 = "R₁"    # Interior radius of the cylindrical shell
    R_2 = "R₂"    # Exterior radius of the cylindrical shell
    x_hat = "x̂"   # Unit vector in the x-direction

    # --- Part 1: Required Permeability (mu) ---
    # From the boundary conditions, we derive the condition (mu/mu_0)^2 = 1.
    # The problem asks to exclude the trivial case mu = mu_0.
    # This leaves the non-trivial solution mu = -mu_0.
    
    print("1. Required Permeability (μ)")
    print("---------------------------------")
    print("For the external field to remain undistorted, the permeability of the shell material must satisfy:")
    # The equation shows the final required relationship between mu and mu_0.
    print(f"μ = -1 * {mu_0}")
    print("\n")

    # --- Part 2: Magnetic Field in the Interior Region (H_int) ---
    # With mu = -mu_0, the magnetic field in the interior region (rho < R1)
    # is found to be uniform and aligned with the external field.
    
    print("2. Magnetic Field in the Interior Region (H_int)")
    print("--------------------------------------------------")
    print("Under this condition, the magnetic field inside the shell is uniform and given by the formula:")
    # The final equation for the interior magnetic field vector.
    # The ^2 indicates squaring of the radii ratio.
    print(f"H_int = {H_0} * ({R_2} / {R_1})^2 * {x_hat}")


# Execute the function to print the solution.
solve_magnetostatics_problem()