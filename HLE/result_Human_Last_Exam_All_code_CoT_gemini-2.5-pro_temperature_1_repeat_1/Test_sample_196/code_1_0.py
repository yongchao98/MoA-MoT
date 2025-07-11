import sympy as sp

def calculate_enclosed_volume():
    """
    Solves for the volume enclosed by a cone S1 and an ellipsoid S2.

    The code follows these steps:
    1. Defines the equations for the ellipsoid and the cone.
    2. Sets up the volume integral using the washer method for volumes of revolution.
    3. Simplifies the integrand.
    4. Calculates the definite integral to find the volume.
    5. Prints a detailed step-by-step evaluation of the final integral calculation.
    """
    y = sp.Symbol('y')
    pi = sp.pi

    # Step 1: Define the radii of the surfaces of revolution
    # The ellipsoid S2 is (x^2)/3 + (y^2)/4 + (z^2)/3 = 1.
    # For a cross-section at height y, the radius r is given by r^2 = x^2 + z^2.
    # So, r^2 / 3 + y^2 / 4 = 1 => r_ellipsoid^2 = 3 * (1 - y^2 / 4)
    r_sq_ellipsoid = 3 * (1 - y**2 / 4)

    # The cone S1 has its vertex at (0, 4, 0) and is tangent to the ellipsoid.
    # The tangency occurs at the plane y=1.
    # The equation of a cone with vertex (0,4,0) and axis along the y-axis is x^2 + z^2 = k*(y-4)^2.
    # At y=1, the radius of the ellipsoid is r^2 = 3*(1 - 1/4) = 9/4.
    # So, 9/4 = k*(1-4)^2 = 9k => k = 1/4.
    # The cone's radius squared is r_cone^2 = (1/4)*(y-4)^2.
    r_sq_cone = (1/4) * (y - 4)**2

    # Step 2: Set up the volume integral
    # The volume is enclosed between the intersection at y=1 and the top of the ellipsoid at y=2.
    # In this interval [1, 2], the cone's radius is larger than the ellipsoid's.
    # V = Integral from 1 to 2 of pi * (r_cone^2 - r_ellipsoid^2) dy
    y_lower, y_upper = 1, 2

    integrand = pi * (r_sq_cone - r_sq_ellipsoid)
    integrand_simplified = sp.simplify(integrand)

    print("--- Volume Calculation ---")
    print(f"The volume is the integral of the difference in cross-sectional areas from y={y_lower} to y={y_upper}.")
    print(f"V = Integral( pi * ( (r_cone^2) - (r_ellipsoid^2) ) ) dy")
    print(f"V = Integral( {integrand} ) dy")
    print(f"Simplifying the integrand, we get:")
    print(f"V = Integral( {integrand_simplified} ) dy")
    print("-" * 30)

    # Step 3: Calculate the definite integral and show steps
    # Antiderivative of pi*(y-1)**2 is pi*(y-1)**3/3
    antiderivative_func = sp.integrate(integrand_simplified, y)
    
    print("The antiderivative of the integrand is:")
    print(f"F(y) = {antiderivative_func}")
    print("\nNow, we evaluate the definite integral from y=1 to y=2:")
    
    # We will manually print the steps for clarity as requested.
    print("\nV = pi * [y^3/3 - y^2 + y] evaluated from 1 to 2")
    print("V = pi * ( ( (2)**3/3 - (2)**2 + 2 ) - ( (1)**3/3 - (1)**2 + 1 ) )")
    print("V = pi * ( ( 8/3 - 4 + 2 ) - ( 1/3 - 1 + 1 ) )")
    print("V = pi * ( ( 8/3 - 2 ) - ( 1/3 ) )")
    print("V = pi * ( ( 8/3 - 6/3 ) - 1/3 )")
    print("V = pi * ( 2/3 - 1/3 )")
    print("V = pi * ( 1/3 )")
    
    # Final result
    final_volume = sp.integrate(integrand_simplified, (y, y_lower, y_upper))
    
    print("-" * 30)
    print(f"The final calculated volume is: {final_volume}")
    print(f"As a decimal, the volume is approximately: {final_volume.evalf()}")

if __name__ == '__main__':
    calculate_enclosed_volume()