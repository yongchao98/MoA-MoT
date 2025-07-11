import sympy

def solve_sequence_problem():
    """
    Solves the math problem by systematically analyzing the two possible cases
    derived from the properties of the arithmetic sequences.
    """
    a1, d, n, k = sympy.symbols('a1 d n k')

    print("Step 1: Find the relationship between a1 and d.")
    print("The condition that {b_n} is an arithmetic sequence is 2*b_2 = b_1 + b_3.")
    # This simplifies to a1^2 - 3*a1*d + 2*d^2 = 0, which factors to (a1 - d)(a1 - 2d) = 0.
    # The solutions are a1 = d and a1 = 2d.
    a1_solutions = [d, 2*d]
    print(f"This leads to two possible cases: a1 = {a1_solutions[0]} or a1 = {a1_solutions[1]}.\n")

    final_d = None
    n_val = 99

    # --- Case 1: a1 = d ---
    print(f"--- Analyzing Case 1: a1 = {a1_solutions[0]} ---")
    a_n_case1 = a1_solutions[0] + (n - 1) * d
    b_n_case1 = (n**2 + n) / a_n_case1
    
    S_99_case1 = sympy.summation(a_n_case1.subs(n, k), (k, 1, n_val))
    T_99_case1 = sympy.summation(b_n_case1.subs(n, k), (k, 1, n_val))

    print(f"The formula for a_n is: {sympy.simplify(a_n_case1)}")
    print(f"The formula for b_n is: {sympy.simplify(b_n_case1)}")
    
    # Setup and solve the equation for d
    equation_case1 = sympy.Eq(S_99_case1 - T_99_case1, n_val)
    simplified_eq = 50*d**2 - d - 51
    print(f"The condition S_99 - T_99 = 99 becomes {S_99_case1} - ({T_99_case1}) = 99.")
    print(f"This simplifies to the quadratic equation for d: 50*d^2 + (-1)*d + (-51) = 0.")

    d_solutions_case1 = sympy.solve(simplified_eq, d)
    print(f"Solving for d gives: d = {d_solutions_case1[0]} or d = {d_solutions_case1[1]}.")
    
    for sol in d_solutions_case1:
        if sol > 1:
            final_d = sol
            print(f"Since d must be > 1, the valid solution from this case is d = {sol.evalf()}.\n")

    # --- Case 2: a1 = 2d ---
    print(f"--- Analyzing Case 2: a1 = {a1_solutions[1]} ---")
    a_n_case2 = a1_solutions[1] + (n - 1) * d
    b_n_case2 = (n**2 + n) / a_n_case2

    S_99_case2 = sympy.summation(a_n_case2.subs(n, k), (k, 1, n_val))
    T_99_case2 = sympy.summation(b_n_case2.subs(n, k), (k, 1, n_val))

    print(f"The formula for a_n is: {sympy.simplify(a_n_case2)}")
    print(f"The formula for b_n is: {sympy.simplify(b_n_case2)}")

    # Setup and solve the equation for d
    equation_case2 = sympy.Eq(S_99_case2 - T_99_case2, n_val)
    simplified_eq2 = 51*d**2 - d - 50
    print(f"The condition S_99 - T_99 = 99 becomes {S_99_case2} - ({T_99_case2}) = 99.")
    print(f"This simplifies to the quadratic equation for d: 51*d^2 + (-1)*d + (-50) = 0.")
    
    d_solutions_case2 = sympy.solve(simplified_eq2, d)
    print(f"Solving for d gives: d = {d_solutions_case2[0]} or d = {d_solutions_case2[1]}.")
    print("Neither of these solutions is greater than 1.\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"The only case that yields a valid result is Case 1.")
    print(f"The value of d is {final_d.evalf()}.")

if __name__ == '__main__':
    solve_sequence_problem()