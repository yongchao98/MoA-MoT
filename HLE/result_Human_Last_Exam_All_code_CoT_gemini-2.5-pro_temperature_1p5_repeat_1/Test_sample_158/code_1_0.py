import sympy

def solve_magnetostatics_problem():
    """
    Solves for the required permeability of a cylindrical shell and the
    resulting interior magnetic field under specific conditions.

    This function outlines the symbolic solution to the problem, explaining
    the physics and mathematical steps that lead to the final answer.
    """
    # Define symbolic variables for the problem
    R1, R2, H0, mu0 = sympy.symbols('R1 R2 H0 mu_0', real=True, positive=True)
    mu = sympy.symbols('mu', real=True) # Permeability of the shell

    # --- Step 1: Theoretical Derivation Summary ---
    # The problem is solved using the magnetic scalar potential V_m, where H = -grad(V_m).
    # The potential must satisfy Laplace's equation in each of the three regions.
    # The key condition is that the external field H = H0*x_hat is undisturbed for rho > R2.
    # This sets the potential in the exterior region to V_ext = -H0 * rho * cos(phi).
    # Applying boundary conditions (continuity of tangential H and normal B) at rho=R1 and rho=R2
    # leads to a system of equations for the potential coefficients and the permeability mu.

    # --- Step 2: The Permeability Condition ---
    # Solving this system of equations leads to the following condition:
    # (mu**2 - mu0**2) * (R2**2 / R1**2) = (mu**2 - mu0**2)
    #
    # This can be rewritten as:
    # (mu**2 - mu0**2) * ( (R2/R1)**2 - 1 ) = 0
    #
    # For a shell of finite thickness, R2 != R1, so the term ((R2/R1)**2 - 1) is non-zero.
    # Therefore, the equation is only satisfied if (mu**2 - mu0**2) = 0.
    # This gives two possible solutions for mu:
    #   1. mu = mu0 (The trivial case, which the problem explicitly excludes)
    #   2. mu = -mu0 (The non-trivial mathematical solution)

    required_mu = -mu0
    relative_permeability = required_mu / mu0

    print("--- Required Permeability of the Shell Material ---")
    print("To prevent distortion of the external magnetic field, the permeability of the shell, mu, must satisfy the relation:")
    print(f"mu = {required_mu}")
    print("The relative permeability, mu_r = mu / mu_0, is therefore:")
    print(f"mu_r = {relative_permeability}")
    print("\nNote: This is a mathematical solution to the boundary value problem. A negative permeability is characteristic of metamaterials at specific frequencies, but typically not of conventional materials in magnetostatics.\n")


    # --- Step 3: Calculation of the Interior Magnetic Field (H_int) ---
    # With mu = -mu0, the coefficients for the potentials can be solved.
    # The potential in the interior region (rho < R1) has the form:
    # V_int = A * rho * cos(phi)
    # Solving the boundary condition equations for A with mu = -mu0 yields:
    # A = -H0 * (R2/R1)**2

    A = -H0 * (R2/R1)**2

    # The magnetic field inside is H_int = -grad(V_int).
    # In Cartesian coordinates, V_int can be written as A*x.
    # So, H_int = -grad(A*x) = -A * x_hat.
    H_int_x = -A

    print("--- Magnetic Field in the Interior Region (rho < R1) ---")
    print("The magnetic field H_int in the interior region is uniform and given by:")
    # We construct a string representation of the final vector expression
    h_int_str = f"H_int = ({sympy.pretty(H_int_x, use_unicode=False)}) * x_hat"
    print(h_int_str)
    print("\nThis means the field inside is also in the x-direction, but its magnitude is magnified by a factor of (R2/R1)^2.")


solve_magnetostatics_problem()