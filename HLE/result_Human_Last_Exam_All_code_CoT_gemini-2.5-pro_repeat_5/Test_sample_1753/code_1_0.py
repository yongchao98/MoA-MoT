import sympy
from sympy import sin, cos, pi, sqrt, acos, Symbol, Eq, solve, pretty_print

def solve_arc_length_problem():
    """
    Solves for the constant 'a' given the parametric arc length problem.
    """
    # Step 1: Define variables and parametric equations
    print("1. Define the parametric equations and variables:")
    t = Symbol('t')
    a = Symbol('a', positive=True, real=True)
    x = cos(t)**3
    y = sin(t)**3
    print(f"   x(t) = {x}")
    print(f"   y(t) = {y}\n")

    # Step 2: Calculate derivatives
    print("2. Calculate the derivatives with respect to t:")
    dx_dt = sympy.diff(x, t)
    dy_dt = sympy.diff(y, t)
    print(f"   dx/dt = {dx_dt}")
    print(f"   dy/dt = {dy_dt}\n")

    # Step 3: Simplify the arc length integrand sqrt((dx/dt)^2 + (dy/dt)^2)
    print("3. Simplify the arc length integrand:")
    integrand_squared = sympy.simplify(dx_dt**2 + dy_dt**2)
    simplified_integrand = sympy.simplify(sqrt(integrand_squared))
    print(f"   (dx/dt)^2 + (dy/dt)^2 = {integrand_squared}")
    print(f"   The simplified integrand is: {simplified_integrand}\n")

    # Step 4: Determine the integration limits from the condition 0 <= x <= a
    print("4. Determine integration limits:")
    print("   The arc is defined by 0 <= x <= a, which means 0 <= cos(t)^3 <= a.")
    print("   This implies 0 <= cos(t) <= a^(1/3).")
    print("   This condition defines two symmetric pieces of the curve corresponding to t in")
    print("   [-pi/2, -arccos(a^(1/3))] and [arccos(a^(1/3)), pi/2].\n")

    # Step 5: Set up and calculate the total arc length integral
    print("5. Calculate the total arc length in terms of 'a':")
    # We calculate the length of the piece in the first quadrant and double it.
    # In the first quadrant, the integrand is 3*sin(t)*cos(t).
    integrand_1st_quad = 3 * sin(t) * cos(t)
    t_a = acos(a**(sympy.S(1)/3))
    
    # Integrate to find the length of one piece
    length_1_piece = sympy.integrate(integrand_1st_quad, (t, t_a, pi/2))
    print(f"   The length of one piece (from t=arccos(a^(1/3)) to t=pi/2) is: {length_1_piece}")
    
    # Total length is twice the length of one piece
    total_length_expr = 2 * length_1_piece
    print(f"   The total length of the two pieces is: 2 * ({length_1_piece}) = {total_length_expr}\n")

    # Step 6: Solve for 'a'
    print("6. Set the total length equal to 3/2 and solve for 'a':")
    given_length = sympy.S(3)/2
    equation = Eq(total_length_expr, given_length)
    
    print("   The equation is:")
    # Print the equation with all the numbers
    pretty_print(equation)

    # Solve the equation
    solution = solve(equation, a)
    final_a = solution[0]
    
    print("\n   Solving the equation step-by-step:")
    print(f"   {total_length_expr} = 3/2")
    print("   a^(2/3) = (3/2) / 3")
    print("   a^(2/3) = 1/2")
    print("   a^2 = (1/2)^3")
    print("   a^2 = 1/8")
    print("   a = sqrt(1/8)")

    print(f"\n   The final value for a is: {final_a}")

if __name__ == '__main__':
    solve_arc_length_problem()

# The final answer is sqrt(2)/4
<<<sqrt(2)/4>>>