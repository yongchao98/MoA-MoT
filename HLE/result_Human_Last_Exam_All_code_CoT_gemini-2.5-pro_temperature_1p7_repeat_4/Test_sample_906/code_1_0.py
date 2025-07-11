import sympy as sp

def solve_steady_state_probability():
    """
    This function finds the steady-state probability pi_0 for the given birth-death process
    using symbolic mathematics.

    The relationship between pi_n and pi_0 is pi_n = pi_0 * (rho^n / n!).
    We use the normalization condition Sum(pi_n for n=0 to infinity) = 1 to solve for pi_0.
    """
    # Define the symbolic variables
    # rho represents the ratio lambda / mu
    rho = sp.Symbol('rho')
    # n is the state index for the summation
    n = sp.Symbol('n', integer=True, nonnegative=True)
    # pi_0 is the steady-state probability of being in state 0, which we want to find
    pi_0_sym = sp.Symbol('pi_0')

    # The term for the sum is pi_n without the pi_0 factor
    # This is rho^n / n!
    series_term = (rho**n) / sp.factorial(n)

    # The normalization condition is Sum(pi_n) = 1, which means
    # pi_0 * Sum(rho^n / n!) = 1.
    # We first calculate the sum of the series.
    series_sum = sp.Sum(series_term, (n, 0, sp.oo)).doit()
    
    # We set up the equation to solve for pi_0: pi_0 * series_sum = 1
    equation = sp.Eq(pi_0_sym * series_sum, 1)

    # Solve the equation for pi_0
    solution = sp.solve(equation, pi_0_sym)

    # The result is a list, so we take the first element.
    pi_0_expression = solution[0]
    
    # Print the final result as a clear equation.
    # The components of the equation are pi_0, e, and -rho.
    final_equation = sp.Eq(pi_0_sym, pi_0_expression)
    print("The derived steady-state probability pi_0 is given by the equation:")
    
    # The pprint function provides a nicely formatted output
    sp.pprint(final_equation, use_unicode=False)

solve_steady_state_probability()