import sympy as sp

def solve_fluid_interface_equation():
    """
    This function derives and prints the coefficients A(r) and B(r)
    for the governing linear equation of the fluid interface.

    The derivation is based on the linearized Young-Laplace equation for
    an axisymmetric fluid interface in cylindrical coordinates.
    """

    # Define symbolic variables for the derivation
    r = sp.Symbol('r', positive=True)  # Radial position
    gamma = sp.Symbol('gamma', positive=True)  # Surface tension
    xi = sp.Function('xi')(r)  # Fluid interface displacement

    print("Step 1: The governing equation for the fluid interface shape xi(r) is based on a pressure balance.")
    print("The balance is between surface tension pressure and electrostatic pressure (encapsulated in C).")
    print("ΔP_tension + C(r, xi) = 0")
    print("-" * 30)

    print("Step 2: The surface tension pressure is given by the Young-Laplace equation, ΔP_tension = gamma * kappa,")
    print("where gamma is the surface tension and kappa is the mean curvature of the interface.")
    print("-" * 30)

    print("Step 3: For an axisymmetric interface xi(r) in cylindrical coordinates, the linearized mean curvature (for small slopes) is:")
    # The linearized curvature kappa is d^2(xi)/dr^2 + (1/r) * d(xi)/dr
    xi_prime = sp.Derivative(xi, r)
    xi_double_prime = sp.Derivative(xi, r, 2)
    kappa_linearized = xi_double_prime + (1/r) * xi_prime
    print(f"kappa ≈ {kappa_linearized}")
    print("-" * 30)

    print("Step 4: The linearized pressure balance equation is therefore:")
    governing_equation = gamma * kappa_linearized + sp.Symbol('C')
    print(f"gamma * ({kappa_linearized}) + C(r, xi) = 0")
    print("Expanding this gives:")
    # We use sp.expand to make it look nicer, though not strictly needed here
    expanded_eq = sp.expand(gamma * kappa_linearized) + sp.Symbol('C')
    print(f"{expanded_eq} + C(r, xi) = 0")
    print("-" * 30)


    print("Step 5: Compare this to the general form provided in the problem:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
    print("-" * 30)

    print("Step 6: By matching the terms, we can identify A(r) and B(r).")
    A_r = gamma
    B_r = gamma / r
    print(f"The coefficient of the second derivative, A(r), is: {A_r}")
    print(f"The coefficient of the first derivative, B(r), is: {B_r}")
    print("-" * 30)
    
    # Final formatted output of the expressions
    print("Final Expressions:")
    print(f"A(r) = {sp.printing.pretty(A_r)}")
    print(f"B(r) = {sp.printing.pretty(B_r)}")


if __name__ == '__main__':
    solve_fluid_interface_equation()
