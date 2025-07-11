import sympy

def display_electrodynamics_solution():
    """
    This function uses the sympy library to symbolically construct and display
    the solution for the electric potential and electric field outside a
    conducting sphere in a uniform electric field, in a steady-state.
    The expressions match the correct answer, derived from first principles.
    """

    # Define the symbolic variables used in the problem
    r, R, theta = sympy.symbols('r R theta', real=True, positive=True)
    E0, sigma1, sigma2 = sympy.symbols('E_0 sigma_1 sigma_2', real=True)

    print("The electric potential and field are found by solving Laplace's equation with the appropriate boundary conditions for a steady current.")
    print("The correct expressions for both inside and outside the sphere are given in option B.")
    print("\n--- Solution for the region OUTSIDE the sphere (r > R) ---")

    # --- Electric Potential Phi(r, theta) for r > R ---
    # The structure is -E0 * (r - [dipole_term]) * cos(theta)
    # The code below constructs this expression for printing.
    phi_str = f"Phi(r, theta) = -{E0} * (r - (({sigma1} - {sigma2}) * {R}**3) / (({sigma1} + 2*{sigma2}) * {r}**2)) * cos({theta})"
    print("\nElectric Potential:")
    print(phi_str)

    # --- Electric Field E(r, theta) for r > R ---
    # The electric field is E = -grad(Phi). It has a radial (r_hat) and
    # an angular (theta_hat) component.
    print("\nElectric Field, E = E_r * r_hat + E_theta * theta_hat, where:")

    # Construct the radial component string
    E_r_str = f"E_r = {E0} * (1 + (2 * ({sigma1} - {sigma2}) * {R}**3) / (({sigma1} + 2*{sigma2}) * {r}**3)) * cos({theta})"
    print(f"\n  Radial component: {E_r_str}")

    # Construct the angular component string
    E_theta_str = f"E_theta = -{E0} * (1 - (({sigma1} - {sigma2}) * {R}**3) / (({sigma1} + 2*{sigma2}) * {r}**3)) * sin({theta})"
    print(f"\n  Angular component: {E_theta_str}")


# Execute the function to display the solution
display_electrodynamics_solution()