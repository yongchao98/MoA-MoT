import sympy as sp

def solve_steady_state_pi0():
    """
    This function symbolically derives the steady-state probability pi_0 for the given
    birth-death process.

    The process has rates:
    - Birth rate (i -> i+1): lambda_i = lambda / (i+1)
    - Death rate (i -> i-1): mu_i = mu

    The function uses the normalization condition sum(pi_k) = 1, where
    pi_k = (rho^k / k!) * pi_0 and rho = lambda / mu.
    """
    # Define the symbols used in the derivation
    # k is the state index (non-negative integer)
    # rho is the traffic intensity (positive real number)
    # pi_0 is the steady-state probability of being in state 0 (positive real number)
    k = sp.symbols('k', integer=True, nonnegative=True)
    rho = sp.symbols('rho', real=True, positive=True)
    pi_0_sym = sp.symbols('pi_0', real=True, positive=True)

    # The relationship pi_k = (rho^k / k!) * pi_0 is derived from the balance equations.
    # We define the term inside the summation for the normalization condition.
    summation_term = rho**k / sp.factorial(k)

    # The normalization condition is: pi_0 * sum_{k=0 to inf} (rho^k / k!) = 1
    # We use sympy to compute the infinite sum, which is the Taylor series for exp(rho).
    total_sum = sp.summation(summation_term, (k, 0, sp.oo))

    # We form the equation: pi_0 * total_sum = 1
    equation = sp.Eq(pi_0_sym * total_sum, 1)

    # We solve this equation for pi_0
    solution = sp.solve(equation, pi_0_sym)

    # The result is a list containing the expression for pi_0
    pi_0_expression = solution[0]

    # Print the final result in a clear format.
    # The final equation is pi_0 = e^(-rho).
    # We explicitly show the numbers involved (0, e, -1) as per the instructions.
    print("The derived steady-state probability for state 0 is:")
    print(f"pi_0 = {pi_0_expression}")
    print("\nWhich can be written with all numbers explicit:")
    final_equation_str = f"pi_0 = e^(-1 * rho)"
    print(final_equation_str)

# Run the function to get the result
solve_steady_state_pi0()