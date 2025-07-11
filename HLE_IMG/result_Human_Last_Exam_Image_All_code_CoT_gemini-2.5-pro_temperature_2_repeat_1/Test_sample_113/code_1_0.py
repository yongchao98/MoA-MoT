import sympy

def solve_fractal_area():
    """
    Calculates the total area of all circles in the fractal pattern as n -> infinity.
    """
    # Step 1: Define initial parameters and calculate the area of the first circle (A0)
    print("--- Step 1: Calculate the area of the first circle (A₀) ---")
    w0 = sympy.Integer(6)
    h0 = sympy.Integer(8)
    
    # The radius of the first circle is one-third of the width
    r0 = w0 / 3
    
    # The area is π * r²
    a0 = sympy.pi * r0**2
    
    print(f"The initial rectangle has width w₀ = {w0} and height h₀ = {h0}.")
    print(f"The radius of the first circle is r₀ = w₀ / 3 = {w0}/3 = {r0}.")
    print(f"The area of the first circle, which is the first term of our series, is A₀ = π * r₀² = {a0}.")
    print("\n")

    # Step 2: Determine the linear scaling factor (s)
    print("--- Step 2: Determine the scaling factor (s) ---")
    # We find the dimensions of the new rectangles formed in the corners.
    # Center the rectangle at (0,0). The diagonal equation is y = (h₀/w₀)x.
    # The circle equation is x² + y² = r₀².
    # We solve for the intersection point's x-coordinate.
    x = sympy.Symbol('x')
    eq = sympy.Eq(x**2 * (1 + (h0/w0)**2), r0**2)
    solutions = sympy.solve(eq, x)
    # We take the positive solution for the x-coordinate in the first quadrant
    x_intersect = max(solutions)

    # The width of the new rectangle in the top-right corner is w₁ = w₀/2 - x_intersect
    w1 = w0/2 - x_intersect
    
    # The linear scaling factor 's' is the ratio of the new width to the original width.
    s = w1 / w0

    print(f"The diagonals and the first circle intersect at x = ±{x_intersect}.")
    print(f"The width of a new smaller rectangle is w₁ = w₀/2 - x_intersect = {w0/2} - {x_intersect} = {w1}.")
    print(f"The linear scaling factor for dimensions is s = w₁ / w₀ = {w1}/{w0} = {s}.")
    print("\n")

    # Step 3: Find the common ratio (r) of the geometric series
    print("--- Step 3: Determine the common ratio (r) of the geometric series ---")
    # At each step, we add 4 new circles. The area of each circle is scaled by s².
    # Thus, the total area added at each step is multiplied by a common ratio r = 4 * s².
    r_geom = 4 * s**2

    print("At each step, 4 new circles are added for each smaller rectangle.")
    print("The area of each of these smaller circles is scaled by s².")
    print(f"The common ratio of the total area added is r = 4 * s² = 4 * ({s})² = {r_geom}.")
    print("\n")

    # Step 4: Sum the infinite geometric series for the total area
    print("--- Step 4: Calculate the total area (S) ---")
    # The total area S is the sum of the infinite geometric series S = A₀ + A₀*r + A₀*r² + ...
    # The formula for the sum is S = A₀ / (1 - r)
    total_area = a0 / (1 - r_geom)

    print("The total area S is the sum of the infinite geometric series S = A₀ / (1 - r).")
    print(f"Plugging in the numbers, the final equation is:")
    print(f"S = ({a0}) / (1 - ({r_geom}))")
    
    one_minus_r = 1 - r_geom
    print(f"S = ({a0}) / ({one_minus_r})")
    print(f"S = {total_area}")
    
    print("\n--- Final Answer ---")
    print(f"The limit of the total area of all circles is {total_area}.")
    print(f"Numerically, this is approximately {total_area.evalf()}.")

solve_fractal_area()