import sympy

def solve_sequence_problem():
    """
    Solves the math problem using symbolic computation with sympy.
    """
    # 1. Define symbolic variables
    a1, d, n = sympy.symbols('a1 d n')
    N = 99

    # 2. Express a_n and b_n
    a_n = a1 + (n - 1) * d
    b_n_expr = (n**2 + n) / a_n

    # 3. Use the arithmetic sequence property for b_n (2*b_2 = b_1 + b_3)
    # to find the relationship between a1 and d.
    # This simplifies to a1^2 - 3*a1*d + 2*d^2 = 0, which factors to (a1-d)(a1-2d)=0.
    # The solutions are a1=d and a1=2d.
    print("From the condition that {b_n} is an arithmetic sequence, we find two possible cases for the first term a1:")
    print("Case 1: a1 = d")
    print("Case 2: a1 = 2*d\n")

    # This variable will store the final valid result for d
    final_d = None

    # --- Case 1: a1 = d ---
    print("--- Analyzing Case 1: a1 = d ---")
    a_n_1 = sympy.simplify(a_n.subs(a1, d))
    b_n_1 = sympy.simplify(b_n_expr.subs(a1, d))
    
    # Calculate sums S_99 and T_99
    S99_1 = sympy.summation(a_n_1, (n, 1, N))
    T99_1 = sympy.summation(b_n_1, (n, 1, N))

    print(f"In this case, a_n = {a_n_1} and b_n = {b_n_1}")
    print(f"The sum S_99 = {S99_1}")
    print(f"The sum T_99 = {T99_1}")
    
    # Formulate and print the equation S_99 - T_99 = 99
    # We want to print the equation with its numerical coefficients.
    s99_coeff = S99_1.coeff(d)
    t99_coeff = sympy.simplify(T99_1 * d).as_numer_denom()[0]
    print(f"The final equation S_99 - T_99 = 99 becomes:")
    print(f"{s99_coeff}*d - {t99_coeff}/d = 99")
    
    # Solve the equation for d
    eq1 = sympy.Eq(S99_1 - T99_1, 99)
    solutions1 = sympy.solve(eq1, d)
    
    print(f"Solving for d gives: {solutions1}")
    # Filter solutions based on d > 1
    for sol in solutions1:
        if sol.is_real and sol > 1:
            final_d = sol
            print(f"Found a valid solution: d = {sol.evalf()}\n")

    # --- Case 2: a1 = 2*d ---
    print("--- Analyzing Case 2: a1 = 2*d ---")
    a_n_2 = sympy.simplify(a_n.subs(a1, 2*d))
    b_n_2 = sympy.simplify(b_n_expr.subs(a1, 2*d))

    # Calculate sums S_99 and T_99
    S99_2 = sympy.summation(a_n_2, (n, 1, N))
    T99_2 = sympy.summation(b_n_2, (n, 1, N))

    print(f"In this case, a_n = {a_n_2} and b_n = {b_n_2}")
    print(f"The sum S_99 = {S99_2}")
    print(f"The sum T_99 = {T99_2}")

    # Formulate and print the equation S_99 - T_99 = 99
    s99_coeff_2 = S99_2.coeff(d)
    t99_coeff_2 = sympy.simplify(T99_2 * d).as_numer_denom()[0]
    print(f"The final equation S_99 - T_99 = 99 becomes:")
    print(f"{s99_coeff_2}*d - {t99_coeff_2}/d = 99")

    # Solve the equation for d
    eq2 = sympy.Eq(S99_2 - T99_2, 99)
    solutions2 = sympy.solve(eq2, d)
    
    print(f"Solving for d gives: {solutions2}")
    print("Neither solution in this case satisfies the condition d > 1.\n")
    
    # --- Conclusion ---
    if final_d is not None:
        print("--- Conclusion ---")
        print(f"The only solution that satisfies all conditions is d = {final_d.evalf()}.")

# Execute the solver function
solve_sequence_problem()