import sympy as sp
from sympy import pi, sqrt

def solve_volume_problem():
    """
    This function calculates the volume enclosed by a cone and an ellipsoid
    by following the outlined plan.
    """
    # Define symbolic variables
    x, y, z = sp.symbols('x y z')

    # Step 1: Define the surfaces
    print("Step 1: Define the surfaces and their equations")
    vertex = (0, 4, 0)
    ellipsoid_eq_str = "x^2/3 + y^2/4 + z^2/3 = 1"
    print(f"S1: A cone with vertex at {vertex}.")
    print(f"S2: An ellipsoid defined by the equation: {ellipsoid_eq_str}.")
    print("-" * 50)

    # Step 2: Find the equation of the cone S1
    # Based on the tangency condition (derived by setting the discriminant of the
    # intersection of a line from the vertex with the ellipsoid to zero),
    # the cone's equation is (y - 4)^2 = 4*(x^2 + z^2).
    cone_eq_xz = (y - 4)**2 / 4
    print("Step 2: Determine the equation of the cone S1")
    print("By enforcing the condition that the cone is tangent to the ellipsoid,")
    print(f"we derive its equation as: x^2 + z^2 = (y - {vertex[1]})^2 / 4")
    print("-" * 50)

    # Step 3: Find the plane of tangency
    print("Step 3: Find the plane of tangency between the cone and the ellipsoid")
    print("We substitute the expression for x^2 + z^2 from the cone's equation into the ellipsoid's equation.")
    # The ellipsoid equation can be written as (x^2 + z^2)/3 + y^2/4 = 1
    # Substituting gives: ((y - 4)^2 / 4) / 3 + y^2/4 = 1
    eq_for_y = cone_eq_xz / 3 + y**2 / 4 - 1
    
    # To solve for y, we simplify the equation:
    # (y^2 - 8y + 16)/12 + 3y^2/12 = 12/12
    # y^2 - 8y + 16 + 3y^2 = 12
    # 4y^2 - 8y + 4 = 0
    # y^2 - 2y + 1 = 0 => (y - 1)^2 = 0
    y_intersect_sol = sp.solve(eq_for_y, y)
    y_intersect = y_intersect_sol[0]
    print(f"Solving the equation '({cone_eq_xz}) / 3 + y^2/4 = 1' for y gives:")
    print("4*y^2 - 8*y + 4 = 0, which simplifies to (y - 1)^2 = 0.")
    print(f"Thus, the surfaces are tangent along the plane where y = {y_intersect}.")
    print("-" * 50)

    # Step 4: Calculate the volume of the cone part (V_cone)
    print("Step 4: Calculate the volume of the cone section (V_cone)")
    h_cone = vertex[1] - y_intersect
    # Find the radius of the cone's base at the tangency plane y=1
    r_sq_cone = cone_eq_xz.subs(y, y_intersect)
    r_cone = sp.sqrt(r_sq_cone)
    # Calculate cone volume using V = 1/3 * pi * r^2 * h
    V_cone_expr = (sp.S(1)/3) * pi * r_sq_cone * h_cone
    
    print(f"This part of the volume is a cone with height h = {vertex[1]} - {y_intersect} = {h_cone}.")
    print(f"Its base radius at y={y_intersect} is r = sqrt(({y_intersect} - {vertex[1]})^2 / 4) = {r_cone}.")
    print("The volume is calculated using the formula V_cone = (1/3) * pi * r^2 * h.")
    print(f"V_cone = (1/3) * pi * ({r_cone})^2 * {h_cone} = {V_cone_expr}")
    print("-" * 50)

    # Step 5: Calculate the volume of the ellipsoid cap (V_ellipsoid_cap)
    print("Step 5: Calculate the volume of the ellipsoid cap (V_cap)")
    # The ellipsoid's y-range is [-2, 2]. We integrate from y=1 to y=2.
    y_top_ellipsoid = 2
    # From the ellipsoid equation, solve for the cross-sectional radius R^2 = x^2 + z^2
    # (x^2 + z^2)/3 = 1 - y^2/4 => x^2 + z^2 = 3 * (1 - y^2/4)
    R_sq_ellipsoid = 3 * (1 - y**2 / 4)
    # Integrate the cross-sectional area pi*R^2 from y=1 to y=2
    V_cap_expr = sp.integrate(pi * R_sq_ellipsoid, (y, y_intersect, y_top_ellipsoid))

    print(f"This part is the volume of the ellipsoid for y from {y_intersect} to {y_top_ellipsoid}.")
    print(f"The cross-sectional area at height y is A(y) = pi * (3 * (1 - y^2/4)).")
    print(f"We integrate A(y) to find the volume: V_cap = integral from {y_intersect} to {y_top_ellipsoid} of A(y) dy.")
    print(f"V_cap = {V_cap_expr}")
    print("-" * 50)
    
    # Step 6: Sum the volumes for the final answer
    print("Step 6: Calculate the total enclosed volume")
    V_total_expr = V_cone_expr + V_cap_expr
    V_total_numeric = V_total_expr.evalf()
    
    print("The total volume is the sum of the cone and cap volumes.")
    print(f"Total Volume = V_cone + V_cap")
    print(f"Total Volume = {V_cone_expr} + {V_cap_expr} = {V_total_expr}")
    print(f"The numerical value of the total volume is approximately {V_total_numeric:.5f}")
    
    return V_total_numeric

if __name__ == '__main__':
    final_volume = solve_volume_problem()
    # The final answer is expected in the format <<<answer>>>
    print(f"\nFinal Answer:")
    print(f"<<<{final_volume:.5f}>>>")
