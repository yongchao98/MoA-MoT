import sympy
from sympy import Sum, factorial, oo, Eq, solve, Function, exp

def solve_steady_state_probability():
    """
    This function symbolically derives the steady-state probability pi_0
    for the given birth-death process.
    """
    # Define symbols
    rho = sympy.Symbol('rho')
    n = sympy.Symbol('n', integer=True, nonnegative=True)
    pi_0 = sympy.Symbol('pi_0')

    # Step 1: Express pi_n in terms of pi_0
    # As derived in the explanation, pi_n = pi_0 * (rho^n / n!)
    pi_n_expr = pi_0 * rho**n / factorial(n)
    print("Step 1: The general steady-state probability pi_n is expressed in terms of pi_0 and rho.")
    # Using sympy's pretty print for better formatting
    sympy.pprint(Eq(Function('pi')(n), pi_n_expr), use_unicode=True)
    print("-" * 40)

    # Step 2: Use the normalization condition Sum(pi_n for n=0 to inf) = 1
    total_prob_sum = Sum(pi_n_expr, (n, 0, oo))
    normalization_eq_symbolic = Eq(total_prob_sum, 1)
    print("Step 2: The normalization equation is set up.")
    sympy.pprint(normalization_eq_symbolic, use_unicode=True)
    print("-" * 40)

    # Step 3: Evaluate the infinite sum. The sum is the Taylor series for exp(rho).
    evaluated_sum = total_prob_sum.doit()
    normalization_eq_evaluated = Eq(evaluated_sum, 1)
    print("Step 3: The infinite sum is evaluated.")
    sympy.pprint(normalization_eq_evaluated, use_unicode=True)
    print("-" * 40)

    # Step 4: Solve the equation for pi_0
    solution = solve(normalization_eq_evaluated, pi_0)
    final_answer = solution[0]

    print("Step 4: The equation is solved for pi_0 to find the final expression.")
    final_equation = Eq(pi_0, final_answer)
    
    # Output the final equation as requested
    print("The final equation is:")
    sympy.pprint(final_equation, use_unicode=True)

if __name__ == '__main__':
    solve_steady_state_probability()