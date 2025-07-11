import sympy
from sympy import pi, sqrt, integrate, symbols

def solve_volume_problem():
    """
    This function calculates the volume of the space enclosed by a cone and an ellipsoid
    by following the plan outlined above.
    """
    # Define symbols for our calculations
    y = symbols('y')

    print("Step 1: Determine the geometry of the intersection.")
    # Based on the derivation in the plan, the intersection is a circle at y=1.
    y_intersect = 1
    # The radius of this circle is r=3/2, so r^2 = 9/4.
    r_intersect_sq = sympy.Rational(9, 4)
    r_intersect = sqrt(r_intersect_sq)
    print(f"The surfaces are tangent along a circle in the plane y = {y_intersect}.")
    print(f"The radius of this circle of tangency is r = {r_intersect}.")
    print("-" * 30)

    print("Step 2: Calculate the volume of the cone segment (V_cone).")
    # The cone segment has its vertex at y=4 and its base at y=1.
    h_cone = 4 - y_intersect
    R_base = r_intersect
    # Volume of a cone: V = (1/3) * pi * r^2 * h
    V_cone = sympy.Rational(1, 3) * pi * (R_base**2) * h_cone
    print(f"The cone segment has height h = 4 - {y_intersect} = {h_cone}.")
    print(f"The base radius is R = {R_base}.")
    print("The volume calculation is:")
    print(f"V_cone = (1/3) * π * ({R_base})^2 * {h_cone} = {V_cone}")
    print("-" * 30)

    print("Step 3: Calculate the volume of the ellipsoid cap (V_cap).")
    # This is the volume of the ellipsoid for y from 1 to its top at y=2.
    y_top_ellipsoid = 2
    # The cross-sectional radius squared r(y)^2 is derived from the ellipsoid equation:
    # r(y)^2 = 3 * (1 - y^2/4)
    r_sq_ellipsoid = 3 * (1 - y**2 / 4)
    integrand = pi * r_sq_ellipsoid
    
    # Perform the integration to find the volume of the cap
    V_cap = integrate(integrand, (y, y_intersect, y_top_ellipsoid))
    print(f"The volume of the cap is the integral of the cross-sectional area A(y) = π * (3*(1 - y^2/4)) from y={y_intersect} to y={y_top_ellipsoid}.")
    print(f"V_cap = ∫[from {y_intersect} to {y_top_ellipsoid}]({integrand}) dy = {V_cap}")
    print("-" * 30)

    print("Step 4: Calculate the final enclosed volume.")
    # The final volume is the difference between the cone segment and the ellipsoid cap.
    V_total = V_cone - V_cap
    print("V_enclosed = V_cone - V_cap")
    print("The final equation for the volume is:")
    print(f"Volume = {V_cone} - {V_cap} = {V_total}")
    
    print("\nIn numerical terms:")
    print(f"Volume ≈ {V_total.evalf()}")

solve_volume_problem()