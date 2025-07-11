import sympy as sp

def solve_fluid_interface_equation():
    """
    Derives the coefficients A(r) and B(r) for the governing equation of a fluid interface
    under the influence of an electric field that modifies surface tension.

    The derivation starts from the Euler-Lagrange equation for the surface energy functional,
    assuming a spatially varying surface tension gamma(r).
    """

    # Define r as a symbol.
    r = sp.Symbol('r')

    # Define gamma and xi as functions of r.
    # gamma(r) is the surface tension, modified by the local electric field.
    # xi(r) is the displacement of the fluid interface.
    gamma = sp.Function('gamma')(r)
    xi = sp.Function('xi')(r)

    # Define the derivatives of xi(r).
    xi_r = sp.diff(xi, r)  # First derivative: d(xi)/dr
    xi_rr = sp.diff(xi_r, r) # Second derivative: d^2(xi)/dr^2

    # The surface energy functional (integrand), ignoring constant factors like 2*pi, for small displacements is:
    # L = gamma(r) * r * (1 + 0.5 * (d(xi)/dr)^2)
    L = gamma * r * (1 + sp.Rational(1, 2) * xi_r**2)

    # According to the Euler-Lagrange equation, the surface tension force is derived from L.
    # The term related to the shape of the interface is -d/dr(dL/d(xi_r)).
    # dL/d(xi_r)
    dL_d_xi_r = sp.diff(L, xi_r)

    # d/dr(dL/d(xi_r))
    surface_tension_force_term = sp.diff(dL_d_xi_r, r)

    # Expand the derivative to identify coefficients of xi_rr and xi_r.
    # The expanded form is f(r)*xi_rr + g(r)*xi_r
    expanded_term = sp.expand(surface_tension_force_term.doit())

    # The equation has the form: A(r)*xi_rr + B(r)*xi_r + C(r, xi) = 0
    # The surface tension force term we derived is -(A(r)*xi_rr + B(r)*xi_r)
    # So we can extract A(r) and B(r) as the coefficients of xi_rr and xi_r.
    
    A_r = expanded_term.coeff(xi_rr)
    B_r = expanded_term.coeff(xi_r)

    # The remaining term, C(r, xi(r)), accounts for other forces like electrostatic pressure not
    # captured in the surface tension modification and is not determined here.

    # Print the derived coefficients A(r) and B(r).
    print("The governing linear equation for the interfacial shape xi(r) is of the form:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0\n")
    print("This derivation assumes the electric field influences the surface tension, making it a function of r, gamma(r).")
    print("Based on the Euler-Lagrange formulation for the interface energy, the coefficients are:\n")
    
    # We output the expressions. As there are no numbers in the problem, we output the symbolic results.
    print(f"A(r) = {A_r}")
    print(f"B(r) = {B_r}")

solve_fluid_interface_equation()