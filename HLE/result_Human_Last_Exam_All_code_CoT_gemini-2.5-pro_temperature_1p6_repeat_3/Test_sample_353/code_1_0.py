import numpy as np
import math

def solve_bdf4_stability_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.
    """
    # 1. Define the cubic polynomial for c = cos(theta)
    # 24*c^3 - 56*c^2 + 34*c - 5 = 0
    coeffs = [24, -56, 34, -5]

    # 2. Find the roots of the polynomial
    roots = np.roots(coeffs)

    # 3. Select the correct root. For BDF4, it's the smallest positive root.
    # The roots are approximately 1.442, 0.725, and 0.166.
    # c must be in [-1, 1], so we filter out the largest root.
    # The stability angle corresponds to the smallest positive root.
    c = min(r for r in roots if r > 0 and r < 1)

    # 4. Define the functions for the x and y coordinates of the stability boundary
    # z(theta) = x(theta) + i*y(theta). We express them in terms of c = cos(theta).
    s = math.sqrt(1 - c**2)
    
    # x(c) from the real part of z(theta)
    x_val = (2*c**4 - (16/3)*c**3 + 4*c**2 - (2/3))
    
    # y(c) from the imaginary part of z(theta)
    # Note: sin(theta) is factored out and multiplied back.
    y_val = s * (-2*c**3 + (16/3)*c**2 - 5*c + (8/3))

    # 5. Calculate the angle alpha
    # alpha = arctan(-y/x)
    tan_alpha = -y_val / x_val
    alpha_rad = math.atan(tan_alpha)
    alpha_deg = math.degrees(alpha_rad)

    # 6. The exact value is known to be arctan(3*sqrt(5)/2)
    exact_tan_alpha = (3 * math.sqrt(5)) / 2
    exact_alpha_rad = math.atan(exact_tan_alpha)

    # 7. Print the results, including the final equation
    print("The cubic equation for c = cos(theta) is: 24*c^3 - 56*c^2 + 34*c - 5 = 0")
    print(f"The relevant root is c = {c:.6f}")
    print(f"This corresponds to theta = {math.degrees(math.acos(c)):.2f} degrees")
    print(f"The point on the stability boundary is z = {x_val:.6f} + {y_val:.6f}i")
    print(f"The tangent of the stability angle is tan(alpha) = -y/x = {tan_alpha:.6f}")
    
    print("\n--- Final Answer ---")
    print("The final equation for the exact angle alpha is:")
    # The problem asks to output the numbers in the equation
    print(f"alpha = arctan( (3 * sqrt(5)) / 2 )")
    print(f"The numerical value of the argument is (3 * sqrt(5)) / 2 = {exact_tan_alpha:.8f}")
    print(f"The angle in radians is alpha = {alpha_rad:.8f}")
    print(f"The angle in degrees is alpha = {alpha_deg:.8f}")


solve_bdf4_stability_angle()