import numpy as np
from scipy import integrate

def solve_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone.
    The cone has base radius 3, height 2, and is centered at the origin.
    """
    # Analytical value for comparison
    analytical_result_val = (108 * np.pi) / 35

    print("We are finding the integral of f(x,y,z) = z^2 * (x^2 + y^2) inside a cone of radius 3 and height 2.")
    print("\nStep 1: Convert to cylindrical coordinates.")
    print("x = r*cos(theta), y = r*sin(theta), z = z")
    print("f(r, theta, z) = z^2 * r^2")
    print("The volume element dV = r dr dz d(theta).")
    print("The integrand is f * r = (z^2 * r^2) * r = z^2 * r^3.")

    print("\nStep 2: Define the bounds of integration.")
    print("For theta: 0 to 2*pi")
    print("For z: 0 to 2")
    print("For r: 0 to 3 - (3/2)*z (the edge of the cone)")
    
    print("\nStep 3: Set up the triple integral.")
    print("Integral = Integral from 0 to 2*pi [ Integral from 0 to 2 [ Integral from 0 to (3 - 1.5*z) [ z^2 * r^3 ] dr ] dz ] d(theta)")

    print("\nStep 4: Solve the integral analytically.")
    # Part 1: Innermost integral (with respect to r)
    print("  1. Integrate z^2 * r^3 with respect to r from 0 to 3 - 1.5*z:")
    print("     z^2 * [r^4 / 4] evaluated from 0 to 3 - 1.5*z = (z^2 / 4) * (3 - 1.5*z)^4")
    
    # Part 2: Middle integral (with respect to z)
    print("  2. Integrate the result with respect to z from 0 to 2:")
    # The analytical result of Integral( (z^2 / 4) * (3 - 1.5*z)^4 dz) is 54/35.
    print("     Integral from 0 to 2 of [(z^2 / 4) * (3 - 1.5*z)^4] dz = 54/35")

    # Part 3: Outermost integral (with respect to theta)
    print("  3. Integrate the result with respect to theta from 0 to 2*pi:")
    print("     Integral from 0 to 2*pi of [54/35] d(theta) = (54/35) * [theta] evaluated from 0 to 2*pi")

    # Final result calculation
    numerator_after_z = 54
    denominator_after_z = 35
    theta_factor = 2
    final_numerator = numerator_after_z * theta_factor
    final_denominator = denominator_after_z

    print("\nFinal Calculation:")
    print(f"Result = ({numerator_after_z} / {denominator_after_z}) * ({theta_factor} * pi) = ({final_numerator}*pi) / {final_denominator}")
    print(f"The numerical value is approximately {analytical_result_val:.4f}")

    # Verification with SciPy tplquad
    # The order of arguments for the integrand is (r, z, theta), matching the integration order.
    # But scipy.integrate.tplquad expects the integrand args as (x,y,z) for (inner, middle, outer)
    integrand = lambda r, z, theta: z**2 * r**3
    
    # We map our integration order dr dz d(theta) to tplquad's (x, y, z)
    # theta is the outer variable 'z' for tplquad
    # z is the middle variable 'y' for tplquad
    # r is the inner variable 'x' for tplquad
    
    r_upper_limit = lambda theta, z: 3 - 1.5 * z
    
    # Note: the arguments to the limit functions are (outer, middle) variables
    # (theta, z) in our case.
    numerical_result, error = integrate.tplquad(
        integrand,
        0, 2 * np.pi,           # theta limits (outer)
        lambda theta: 0,        # z lower limit (middle)
        lambda theta: 2,        # z upper limit (middle)
        lambda theta, z: 0,     # r lower limit (inner)
        r_upper_limit           # r upper limit (inner)
    )

    print(f"\nVerification using scipy.integrate.tplquad: {numerical_result:.4f}")

solve_integral()
<<<9.6945>>>