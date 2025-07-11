import sympy
from sympy import pi, integrate, Symbol, simplify

def solve_volume_problem():
    """
    Calculates the volume of the space enclosed by a cone and an ellipsoid
    and prints the step-by-step solution.
    """
    # Define the variable of integration
    y = Symbol('y')

    # Define the bounds of integration, based on the plane of tangency y=1
    # and the top of the ellipsoid y=2.
    y_min, y_max = 1, 2

    # --- Define area functions ---
    # The cross-sectional area of the cone S1 at height y.
    # The cone's equation is x^2 + z^2 = (9/4) * ((4-y)/3)^2 = (4-y)^2 / 4.
    A_cone = pi * (4 - y)**2 / 4

    # The cross-sectional area of the ellipsoid S2 at height y.
    # The ellipsoid's equation gives x^2 + z^2 = 3 * (1 - y^2/4).
    A_ellipsoid = 3 * pi * (1 - y**2 / 4)

    # --- Print problem setup and explanation ---
    print("This script calculates the volume of the space enclosed by a cone S1 and an ellipsoid S2.")
    print("-" * 50)
    print("Step 1: Determine the integration limits.")
    print("The cone S1 with vertex (0,4,0) is tangent to the ellipsoid S2 on the plane y = 1.")
    print("The ellipsoid S2 extends up to y = 2.")
    print(f"Thus, the volume is found by integrating from y = {y_min} to y = {y_max}.")
    print("-" * 50)

    # --- Print integral setup ---
    print("Step 2: Define the volume as an integral of the difference in cross-sectional areas.")
    print("V = Integral from y=1 to y=2 of [A_cone(y) - A_ellipsoid(y)] dy")
    print(f"A_cone(y) = {A_cone}")
    print(f"A_ellipsoid(y) = {A_ellipsoid}")

    # --- Simplify and calculate the integral ---
    integrand = A_cone - A_ellipsoid
    simplified_integrand = simplify(integrand)

    print("\nThe integrand simplifies as follows:")
    print(f"{A_cone} - ({A_ellipsoid})")
    print(f"= {simplified_integrand}")
    print("-" * 50)

    print("Step 3: Calculate the definite integral.")
    antiderivative = integrate(simplified_integrand, y)
    print(f"V = Integral from {y_min} to {y_max} of ({simplified_integrand}) dy")
    print(f"The antiderivative is: [ {antiderivative} ]")
    
    # --- Show the final equation with all numbers ---
    print("\nEvaluating the antiderivative at the limits y=2 and y=1:")
    print(f"V = [ {pi} * ({y_max} - 1)**3 / 3 ] - [ {pi} * ({y_min} - 1)**3 / 3 ]")
    print(f"V = [ {pi} * (1)**3 / 3 ] - [ {pi} * (0)**3 / 3 ]")
    print(f"V = [ {pi}/3 ] - [ 0 ]")
    
    final_volume = integrate(simplified_integrand, (y, y_min, y_max))

    print("-" * 50)
    print(f"The final volume is: {final_volume}")
    print(f"The numerical value is approximately: {final_volume.evalf()}")

solve_volume_problem()