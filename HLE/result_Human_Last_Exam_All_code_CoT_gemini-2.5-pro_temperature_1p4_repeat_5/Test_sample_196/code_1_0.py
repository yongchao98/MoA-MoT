import sympy

def solve_volume_problem():
    """
    This function calculates the volume of the space enclosed by the cone S1 and the ellipsoid S2.
    """
    # Define symbols for calculation
    y = sympy.Symbol('y')
    pi = sympy.pi

    print("Step 1: Define the geometry and find the plane of tangency.")
    # The ellipsoid is x^2/3 + y^2/4 + z^2/3 = 1.
    # The cone has its vertex at (0, 4, 0).
    # The tangent plane at (x0, y0, z0) is x*x0/3 + y*y0/4 + z*z0/3 = 1.
    # Substitute the vertex (0, 4, 0) to find the plane of tangency:
    # 0*x0/3 + 4*y0/4 + 0*z0/3 = 1  =>  y0 = 1.
    y_tangent = 1
    print(f"The plane of tangency is at y = {y_tangent}.\n")

    print("Step 2: Find the radius of the circle of tangency.")
    # Substitute y=1 into the ellipsoid equation: x^2/3 + 1^2/4 + z^2/3 = 1
    # x^2 + z^2 = 3 * (1 - 1/4) = 9/4.
    r_squared = sympy.Rational(9, 4)
    r = sympy.sqrt(r_squared)
    print(f"The circle of tangency has radius squared r^2 = {r_squared}, so radius r = {r}.\n")

    print("Step 3: Calculate the volume of the cone segment (V_cone).")
    # The cone segment has its vertex at y=4 and its base at y=1.
    h_cone = 4 - y_tangent
    # The volume of a cone is (1/3) * pi * r^2 * h.
    V_cone = (sympy.Rational(1, 3)) * pi * r_squared * h_cone
    print(f"The cone has height h = 4 - 1 = {h_cone} and base radius r = {r}.")
    print(f"V_cone = (1/3) * pi * ({r})^2 * {h_cone} = {V_cone}\n")

    print("Step 4: Calculate the volume of the ellipsoid cap (V_cap).")
    # The volume is the integral of the cross-sectional area A(y) from y=1 to y=2.
    # From the ellipsoid equation, the radius squared of a cross-section at height y is r(y)^2 = 3 * (1 - y^2/4).
    # A(y) = pi * r(y)^2.
    cross_sectional_area = pi * 3 * (1 - y**2 / 4)
    # The ellipsoid top is at y=2 (from y^2/4=1).
    y_top_ellipsoid = 2
    V_cap = sympy.integrate(cross_sectional_area, (y, y_tangent, y_top_ellipsoid))
    print(f"The ellipsoid cap is the volume from y={y_tangent} to y={y_top_ellipsoid}.")
    print(f"V_cap = integral from {y_tangent} to {y_top_ellipsoid} of {cross_sectional_area} dy = {V_cap}\n")

    print("Step 5: Calculate the final enclosed volume.")
    # The total enclosed volume is V = V_cone - V_cap.
    total_volume = V_cone - V_cap
    print(f"The final enclosed volume is the difference between the cone segment and the ellipsoid cap.")
    # The final equation as required
    print("Final Equation:")
    print(f"Volume = V_cone - V_cap")
    print(f"Volume = {V_cone} - {V_cap}")
    print(f"Volume = {total_volume}")
    print("\nThe numerical value of the volume is:")
    print(total_volume.evalf())

solve_volume_problem()