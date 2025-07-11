import sympy

def solve_for_t_bounds():
    """
    This function solves for the bounds of t based on the derived inequalities.
    """
    # Define the symbolic variable t
    t = sympy.Symbol('t')

    # The problem requires t to be in the interval [-1, 0)
    domain = sympy.Interval.Lopen(-1, 0)

    # The condition [1/(2t), -1/2] is a subset of [-2, 2t] leads to two inequalities.
    # The numbers in these equations are -2, 1, 2 and -1/2, 2.
    
    # Inequality 1: -2 <= 1/(2*t)
    ineq1 = -2 <= 1 / (2 * t)
    
    # Inequality 2: -1/2 <= 2*t
    ineq2 = sympy.Rational(-1, 2) <= 2 * t

    # Solve each inequality within the valid domain for t
    sol1 = sympy.solveset(ineq1, t, domain=domain)
    sol2 = sympy.solveset(ineq2, t, domain=domain)

    # The final solution for t must satisfy both inequalities, so we find the intersection
    solution_set = sympy.Intersection(sol1, sol2)

    # The problem asks for the lower and upper bounds of the set of valid t.
    if solution_set.is_empty:
        print("No solution found for t.")
        return

    lower_bound = solution_set.inf
    upper_bound = solution_set.sup

    print(f"The first inequality is (-2 <= 1/(2*t)), which for t in [-1, 0) solves to: {sol1}")
    print(f"The second inequality is (-1/2 <= 2*t), which for t in [-1, 0) solves to: {sol2}")
    print(f"The intersection of these solutions is: {solution_set}")
    print(f"Thus, the only valid value for t is {lower_bound}.")
    
    # Print the final answer in the required format
    print(f"{float(lower_bound)} {float(upper_bound)}")

solve_for_t_bounds()