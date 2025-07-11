import sympy

def solve_volume_problem():
    """
    Calculates the volume of the space enclosed by a cone and an ellipsoid.
    """
    # Define symbolic variable y
    y = sympy.Symbol('y')

    # --- Step 1 & 2: Analyze geometry and find tangency ---
    # The ellipsoid is S2: x^2/3 + y^2/4 + z^2/3 = 1. Top is at y=2.
    # The cone S1 has its vertex at (0, 4, 0).
    # The tangency between S1 and S2 occurs on the plane y = 1.
    # At y=1, the intersection is a circle x^2 + z^2 = 9/4, with radius 3/2.

    # --- Step 3: Set up expressions for the radii ---
    # The enclosed volume is bounded by y=1 and y=2 (top of the ellipsoid).
    # We will integrate the area of annular cross-sections from y=1 to y=2.

    # Squared radius of the ellipsoid's cross-section at height y:
    # From S2, let r^2 = x^2 + z^2. Then r^2/3 + y^2/4 = 1.
    r_ellipsoid_sq = 3 * (1 - y**2 / 4)

    # Squared radius of the cone's cross-section at height y:
    # The cone's profile is a line in the r-y plane through (r=0, y=4) and (r=3/2, y=1).
    # The equation is r = (4 - y) / 2.
    r_cone_sq = ((4 - y) / 2)**2

    print("The volume is found by integrating the area of annular cross-sections from y=1 to y=2.")
    print("The area at height y is A(y) = pi * (r_cone^2 - r_ellipsoid^2).\n")

    print("Squared radius of the cone:")
    print(f"r_C^2 = ((4 - y)/2)^2 = {sympy.expand(r_cone_sq)}")
    print("Squared radius of the ellipsoid:")
    print(f"r_E^2 = 3 * (1 - y^2/4) = {sympy.expand(r_ellipsoid_sq)}\n")

    # --- Step 4 & 5: Calculate the integral ---
    # The integrand is the area of the annulus, A(y)
    A_y = sympy.pi * (r_cone_sq - r_ellipsoid_sq)
    integrand_simplified = sympy.simplify(A_y)
    
    print("The area to be integrated is:")
    print(f"A(y) = {integrand_simplified}\n")
    
    # The definite integral for the volume V
    volume_integral = sympy.Integral(integrand_simplified, (y, 1, 2))
    
    print("The volume calculation is the definite integral:")
    print(f"Volume V = {volume_integral}\n")

    # Solve the integral
    antiderivative = sympy.integrate(integrand_simplified, y)
    final_volume = sympy.integrate(integrand_simplified, (y, 1, 2))
    
    print("--- Final Calculation Breakdown ---")
    print(f"The indefinite integral of {integrand_simplified} with respect to y is: {antiderivative}")
    print("We evaluate this from y=1 to y=2:")

    # Show the substitution with numbers
    antiderivative_func = sympy.lambdify(y, antiderivative/sympy.pi)
    term2 = antiderivative_func(2)
    term1 = antiderivative_func(1)
    
    print(f"Volume = \u03c0 * [ (({2})^3/3 - ({2})^2 + ({2})) - (({1})^3/3 - ({1})^2 + ({1})) ]")
    print(f"Volume = \u03c0 * [ ({sympy.Rational(8,3)} - {4} + {2}) - ({sympy.Rational(1,3)} - {1} + {1}) ]")
    print(f"Volume = \u03c0 * [ ({sympy.Rational(8,3) - 2}) - ({sympy.Rational(1,3)}) ]")
    print(f"Volume = \u03c0 * [ ({sympy.Rational(2,3)}) - ({sympy.Rational(1,3)}) ]")
    print(f"Volume = \u03c0 * ({sympy.Rational(1,3)})")
    
    print("\n--- Final Answer ---")
    print(f"The volume of the space is: {final_volume}")

solve_volume_problem()