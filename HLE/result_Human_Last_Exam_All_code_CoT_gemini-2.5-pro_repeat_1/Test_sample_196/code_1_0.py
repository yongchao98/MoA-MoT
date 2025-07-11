import sympy
from sympy import pi

def solve_volume_problem():
    """
    Calculates the volume of the space enclosed by a cone and an ellipsoid.
    """
    # Define symbols
    y = sympy.Symbol('y')
    
    # Parameters derived from the problem description
    # The cone S1 has its vertex at (0,4,0).
    # The ellipsoid S2 is (x^2/3) + (y^2/4) + (z^2/3) = 1.
    # The cone is tangent to the ellipsoid. The plane of tangency for a cone with a vertex
    # at (0,4,0) to this ellipsoid is y=1.
    # The intersection curve is found by substituting y=1 into the ellipsoid's equation:
    # (x^2/3) + (1/4) + (z^2/3) = 1  => x^2 + z^2 = 9/4.
    # This is a circle of radius r = 3/2.
    
    # Cone parameters
    r_cone_base = sympy.Rational(3, 2)
    h_cone = sympy.Integer(3) # Height from vertex y=4 to base y=1

    # Print introduction and setup
    print("Step-by-step calculation of the volume:")
    print("1. Identify the surfaces and their intersection.")
    print("   - Ellipsoid S2: (x^2/3) + (y^2/4) + (z^2/3) = 1")
    print("   - Cone S1: Vertex at (0,4,0), tangent to S2.")
    print("   - The plane of tangency is found to be y = 1.")
    print(f"   - The curve of tangency is the circle x^2 + z^2 = ({r_cone_base})^2 on the plane y=1.")
    
    print("\n2. Define the enclosed volume.")
    print("   The volume is the space between the cone (upper surface) and the ellipsoid (lower surface).")
    print("   It can be calculated as the volume of a cone section minus the volume of an ellipsoid segment.")
    print("   V = V_cone - V_ellipsoid_segment")

    # Step 1: Calculate the volume of the cone part (V_cone)
    # This is a cone with its vertex at y=4 and its base on the plane y=1.
    # Formula for cone volume: V_cone = (1/3) * pi * r^2 * h
    V_cone = sympy.Rational(1, 3) * pi * r_cone_base**2 * h_cone
    
    print("\n3. Calculate the volume of the cone (V_cone).")
    print(f"   The cone has its vertex at y=4 and its base on the plane y=1.")
    print(f"   Height h = 4 - 1 = {h_cone}")
    print(f"   Base radius r = {r_cone_base}")
    print(f"   V_cone = (1/3) * pi * r^2 * h = (1/3) * pi * ({r_cone_base})^2 * {h_cone}")
    print(f"   V_cone = {V_cone}")

    # Step 2: Calculate the volume of the ellipsoid segment (V_ellipsoid_segment)
    # This is the volume of the ellipsoid cap for y >= 1.
    # We integrate the area of circular cross-sections from y=1 to y=2 (top of ellipsoid).
    # From the ellipsoid equation, the radius squared of a cross-section at height y is R(y)^2 = 3 * (1 - y^2/4).
    # The area of the cross-section is A(y) = pi * R(y)^2.
    A_y_expr = 3 * pi * (1 - y**2 / 4)
    V_ellipsoid_segment = sympy.integrate(A_y_expr, (y, 1, 2))
    
    print("\n4. Calculate the volume of the ellipsoid segment (V_ellipsoid_segment).")
    print("   This is the volume of the ellipsoid for y from 1 to 2.")
    print("   We integrate the cross-sectional area A(y) from y=1 to y=2.")
    print(f"   The cross-sectional area is A(y) = {A_y_expr}.")
    print(f"   V_ellipsoid_segment = integral from 1 to 2 of ({A_y_expr}) dy")
    
    # Show integration details
    integral_result = sympy.integrate(A_y_expr, y)
    antiderivative_at_2 = integral_result.subs(y, 2)
    antiderivative_at_1 = integral_result.subs(y, 1)
    print(f"   The integral evaluates to [{integral_result}] from 1 to 2.")
    print(f"   = ({antiderivative_at_2}) - ({antiderivative_at_1})")
    print(f"   = {V_ellipsoid_segment}")

    # Step 3: Calculate the final volume
    V_final = V_cone - V_ellipsoid_segment
    
    print("\n5. Calculate the final volume.")
    print(f"   V = V_cone - V_ellipsoid_segment")
    print(f"   V = {V_cone} - {V_ellipsoid_segment}")
    print(f"   V = {V_final}")

if __name__ == '__main__':
    solve_volume_problem()