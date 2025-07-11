import sympy

def solve_sequence_problem():
    """
    Solves the math problem by analyzing the two possible cases for the sequence a_n.
    """
    d = sympy.Symbol('d')
    final_d_solutions = []

    print("The problem is solved by considering two cases derived from the properties of the sequences.")
    print("The condition that {b_n} is an arithmetic progression leads to (a_1 - d)(a_1 - 2d) = 0.")
    print("This gives two cases: a_1 = d or a_1 = 2d.\n")

    # --- Case 1: a_1 = d ---
    print("--- Analyzing Case 1: a_1 = d ---")
    # In this case, a_n = d + (n-1)d = n*d.
    # b_n = n(n+1) / (n*d) = (n+1)/d.
    # S_99 = sum_{k=1 to 99} k*d = d * 99*100/2 = 4950*d
    # T_99 = sum_{k=1 to 99} (k+1)/d = (1/d) * sum_{j=2 to 100} j = (1/d) * (100*101/2 - 1) = 5049/d
    # The equation S_99 - T_99 = 99 becomes 4950*d - 5049/d = 99.
    # Multiplying by d and rearranging gives: 4950*d^2 - 99*d - 5049 = 0.
    # Dividing by 99 simplifies the equation to: 50*d^2 - d - 51 = 0.
    
    # Define the coefficients of the quadratic equation 50d^2 - d - 51 = 0
    c2, c1, c0 = 50, -1, -51
    print("The equation S_99 - T_99 = 99 leads to the following quadratic equation for d:")
    print(f"{c2} * d^2 + ({c1}) * d + ({c0}) = 0")
    
    # Solve for d
    equation1 = c2*d**2 + c1*d + c0
    d_solutions1 = sympy.solve(equation1, d)
    print(f"Solutions for d in this case are: {d_solutions1}")
    
    # Check the condition d > 1
    for sol in d_solutions1:
        if sol.is_real and sol > 1:
            final_d_solutions.append(sol)
            print(f"Valid solution found: d = {float(sol):.2f}\n")
        else:
            print(f"Solution d = {sol} is rejected (condition: d > 1).\n")


    # --- Case 2: a_1 = 2d ---
    print("--- Analyzing Case 2: a_1 = 2d ---")
    # In this case, a_n = 2d + (n-1)d = (n+1)*d.
    # b_n = n(n+1) / ((n+1)*d) = n/d.
    # S_99 = sum_{k=1 to 99} (k+1)*d = d * sum_{j=2 to 100} j = 5049*d
    # T_99 = sum_{k=1 to 99} k/d = (1/d) * 99*100/2 = 4950/d
    # The equation S_99 - T_99 = 99 becomes 5049*d - 4950/d = 99.
    # Multiplying by d and rearranging gives: 5049*d^2 - 99*d - 4950 = 0.
    # Dividing by 99 simplifies the equation to: 51*d^2 - d - 50 = 0.

    # Define the coefficients of the quadratic equation 51d^2 - d - 50 = 0
    c2, c1, c0 = 51, -1, -50
    print("The equation S_99 - T_99 = 99 leads to the following quadratic equation for d:")
    print(f"{c2} * d^2 + ({c1}) * d + ({c0}) = 0")
    
    # Solve for d
    equation2 = c2*d**2 + c1*d + c0
    d_solutions2 = sympy.solve(equation2, d)
    print(f"Solutions for d in this case are: {d_solutions2}")
    
    # Check the condition d > 1
    for sol in d_solutions2:
        if sol.is_real and sol > 1:
            final_d_solutions.append(sol)
            print(f"Valid solution found: d = {sol}\n")
        else:
            print(f"Solution d = {sol} is rejected (condition: d > 1).\n")

    # --- Conclusion ---
    print("--- Final Result ---")
    if len(final_d_solutions) == 1:
        final_answer = final_d_solutions[0]
        print(f"The only solution that satisfies the condition d > 1 is d = {final_answer}.")
        print(f"The final value of d is {float(final_answer):.2f}")
    elif len(final_d_solutions) == 0:
        print("No solution found that satisfies the condition d > 1.")
    else:
        print("Multiple valid solutions found.")

solve_sequence_problem()