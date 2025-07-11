import math

def solve_hipparchus_problem():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model
    approximating motion around a square.

    The motion of an object around a square centered at the origin can be
    represented by a complex function z(t). A deferent-epicycle model is an
    approximation of z(t) using the two largest terms of its Fourier series:
    z(t) ≈ c_1 * exp(i * ω_0 * t) + c_k * exp(i * k * ω_0 * t)

    Due to the 4-fold symmetry of the square, the Fourier series only contains
    terms with indices n = 1, -3, 5, -7, ... The deferent corresponds to the
    largest term (n=1), and the epicycle to the second-largest (n=-3).

    R is the ratio of the radii of the deferent and epicycle: R = |c_1| / |c_{-3}|.
    phi is the ratio of their frequencies: phi = (-3 * ω_0) / (1 * ω_0) = -3.

    The analytical calculation for R yields the formula: R = 9 * (π + 2) / (3π - 2).
    This script calculates the numerical value.
    """

    # Use the symbolic formula derived from the Fourier analysis.
    pi = math.pi
    R = 9 * (pi + 2) / (3 * pi - 2)

    # The frequency ratio phi is -3 from the indices of the Fourier terms.
    phi = -3.0

    # Output the components of the final answer as requested.
    print(f"The radius ratio is R = 9 * (pi + 2) / (3*pi - 2)")
    print(f"Calculated value of R: {R}")
    print(f"The frequency ratio is phi = -3")
    print(f"Final equation for the ordered pair:")
    print(f"(R, phi) = ({R}, {phi})")

solve_hipparchus_problem()