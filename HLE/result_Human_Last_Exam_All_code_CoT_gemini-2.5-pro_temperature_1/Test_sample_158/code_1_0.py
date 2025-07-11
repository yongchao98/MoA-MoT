import sympy as sp

def solve_magnetostatics_problem():
    """
    This script presents the solution to the cylindrical shell magnetostatics problem.

    It calculates and displays the required permeability of the shell and the
    resulting magnetic field inside the shell, based on the derivation that
    the external magnetic field remains undisturbed.
    """

    # Define symbolic variables for clarity, although the result is numerical.
    mu_0, H_0, R_1, R_2 = sp.symbols('mu_0 H_0 R1 R2')
    x_hat = sp.Symbol('x_hat')

    # --- Part 1: Determine the required permeability ---
    # From the derivation, the relative permeability must satisfy mu_r^2 = 1.
    # Excluding the trivial case mu_r = 1, the solution is mu_r = -1.
    mu_r_val = -1
    
    print("--- Problem Solution ---")
    print("\nPart 1: Required Permeability of the Shell")
    print("To ensure the external magnetic field is not distorted, the relative permeability (mu_r = mu/mu_0) of the shell must be:")
    print(f"mu_r = {mu_r_val}")
    print("Therefore, the permeability of the shell is:")
    print(f"mu = {mu_r_val} * mu_0")

    # --- Part 2: Determine the interior magnetic field ---
    # The derivation for the interior field H_int = H_int_mag * x_hat gives
    # H_int_mag = H0 * (R2/R1)^2.
    # The coefficient is (R2/R1) raised to the power of 2.
    exponent = 2

    print("\nPart 2: Magnetic Field in the Interior Region")
    print("With this permeability, the magnetic field in the interior region (rho < R1) is uniform and given by:")
    # Using an f-string to clearly show the equation structure and the numerical exponent
    print(f"H_int = H0 * (R2/R1)**{exponent} * x_hat")
    print("\nThis result shows that the shell acts as a magnetic field concentrator,")
    print("increasing the field strength in the core without disturbing the field outside.")

# Execute the function to print the solution.
solve_magnetostatics_problem()