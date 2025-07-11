import sympy

def find_common_difference_d():
    """
    This function solves the problem by following the logical steps derived above.
    It uses the sympy library for symbolic mathematics to avoid manual calculation errors.
    """
    # 1. Define symbolic variables
    a1, d, n = sympy.symbols('a1 d n')

    # 2. Set up the expressions for a_n and b_n
    an_expr = a1 + (n - 1) * d
    bn_expr = (n**2 + n) / an_expr

    # 3. Use the arithmetic property of b_n (2*b2 = b1 + b3) to find relation between a1 and d
    b1 = bn_expr.subs(n, 1)
    b2 = bn_expr.subs(n, 2)
    b3 = bn_expr.subs(n, 3)
    
    # Equation derived from 2*b2 = b1 + b3 is a1^2 - 3*a1*d + 2*d^2 = 0
    relation_eq = a1**2 - 3*a1*d + 2*d**2
    a1_solutions = sympy.solve(relation_eq, a1)

    print(f"From the condition that {{b_n}} is an arithmetic sequence, we get two possible cases for a1: {a1_solutions}\n")

    final_d_solution = None

    # 4. Evaluate each case
    for i, a1_sol in enumerate(a1_solutions):
        print(f"--- Analyzing Case {i+1}: a1 = {a1_sol} ---")
        an_case = an_expr.subs(a1, a1_sol).simplify()
        bn_case = sympy.simplify((n**2 + n) / an_case)
        
        # Calculate sums S_99 and T_99 using arithmetic series sum formula: N/2 * (first + last)
        N = 99
        S99 = (N / 2) * (an_case.subs(n, 1) + an_case.subs(n, N))
        T99 = (N / 2) * (bn_case.subs(n, 1) + bn_case.subs(n, N))
        
        # 5. Set up and solve the equation S_99 - T_99 = 99
        final_eq = sympy.Eq(S99 - T99, 99)
        
        # To show the equation clearly, let's form the polynomial
        # (S99 - T99 - 99)*d gives a polynomial in d if the original had a 1/d term.
        poly_expr = sympy.fraction(sympy.cancel(S99 - T99 - 99))[0]
        poly = sympy.Poly(poly_expr, d)
        
        A = poly.coeffs()[0]
        B = poly.coeffs()[1]
        C = poly.coeffs()[2]
        gcd = sympy.gcd(list(poly.coeffs()))
        
        print(f"The equation S_99 - T_99 = 99 simplifies to the quadratic equation:")
        print(f"Equation: {int(A/gcd)}*d**2 + {int(B/gcd)}*d + {int(C/gcd)} = 0")
        print(f"The numbers in the equation are: a={int(A/gcd)}, b={int(B/gcd)}, c={int(C/gcd)}")

        d_sols = sympy.solve(final_eq, d)
        print(f"Solutions for d: {d_sols}")

        # 6. Check for the condition d > 1
        for sol in d_sols:
            if sol.is_real and sol > 1:
                final_d_solution = sol
                print(f"Valid solution found: d = {sol.evalf()}, which is greater than 1.\n")
            else:
                print(f"Solution d = {sol.evalf()} is discarded because it does not satisfy d > 1.\n")
                
    if final_d_solution is not None:
        print("Final Answer:")
        print(final_d_solution.evalf())

find_common_difference_d()