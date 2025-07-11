import sympy as sp

def find_pi_zero():
    """
    Symbolically derives the steady-state probability pi_0 for a birth-death process
    where birth rate is lambda/(i+1) and death rate is mu.

    The final equation derived is pi_0 = e**(-rho), where rho = lambda / mu.
    """
    # Define the symbolic variables
    rho = sp.symbols('rho', positive=True)
    k = sp.symbols('k', integer=True, nonnegative=True)
    pi_0_sym = sp.Symbol('pi_0')

    # From the detailed balance equations, we derive pi_k in terms of pi_0
    # pi_k = pi_0 * (rho**k / k!)
    pi_k = pi_0_sym * (rho**k / sp.factorial(k))

    # The sum of all probabilities must be 1 (Normalization Condition)
    # Sum(pi_k for k=0 to infinity) = 1
    # This leads to: pi_0 * Sum(rho**k / k!) = 1
    sum_term = sp.Sum(rho**k / sp.factorial(k), (k, 0, sp.oo))

    # The sum is the Taylor series for e**rho
    equation = sp.Eq(pi_0_sym * sum_term.doit(), 1)

    # Solve for pi_0
    solution = sp.solve(equation, pi_0_sym)
    pi_0_expression = solution[0]

    # Create the final equation for printing
    final_equation = sp.Eq(sp.Symbol('pi_0'), pi_0_expression)

    # Print the final equation
    # Using sympy's pretty print for a clear mathematical representation
    print("The final equation for the steady-state probability pi_0 is:")
    print(sp.pretty(final_equation, use_unicode=True))


if __name__ == "__main__":
    find_pi_zero()
