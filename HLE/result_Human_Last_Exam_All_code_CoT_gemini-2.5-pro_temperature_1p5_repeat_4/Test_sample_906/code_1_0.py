import sympy
from sympy import symbols, Eq, Sum, oo, factorial, pprint

def find_pi_0():
    """
    This function analytically derives the steady-state probability pi_0
    for a given birth-death process and prints the result.
    
    The rates are:
    - Birth rate at state i: lambda_i = lambda / (i + 1)
    - Death rate at state i: mu_i = mu
    """

    print("Step 1: Define the relationship between steady-state probabilities.")
    print("For this birth-death process, the relationship between pi_n and pi_0 is:")
    print("pi_n = pi_0 * (rho**n / n!), where rho = lambda / mu.\n")

    # Define symbolic variables for the equation
    pi_0 = symbols('pi_0')
    rho = symbols('rho', positive=True)
    n = symbols('n', integer=True, nonnegative=True)

    # Expression for pi_n
    pi_n_expr = pi_0 * rho**n / factorial(n)

    print("Step 2: Apply the normalization condition.")
    print("The sum of all steady-state probabilities must equal 1:")
    print("Sum_{n=0 to inf} pi_n = 1\n")
    
    # Create the equation Sum(pi_n) = 1
    normalization_eq = Eq(Sum(pi_n_expr, (n, 0, oo)), 1)
    
    print("Substituting the expression for pi_n, we get the equation:")
    # Use use_unicode=False for compatibility with all terminals
    pprint(normalization_eq, use_unicode=False)
    print()

    print("Step 3: Solve the equation for pi_0.")
    print("The summation Sum_{n=0 to inf} (rho**n / n!) is the Taylor series for e**rho.")

    # The `doit()` method evaluates the sum
    solved_sum_eq = Eq(pi_0 * Sum(rho**n / factorial(n), (n, 0, oo)).doit(), 1)
    
    print("So, the equation simplifies to:")
    pprint(solved_sum_eq, use_unicode=False)
    print()
    
    # Solve for pi_0
    final_solution = sympy.solve(solved_sum_eq, pi_0)[0]

    print("Step 4: The final expression for pi_0 is:")
    final_equation = Eq(symbols('pi_0'), final_solution)

    # Pretty print the final equation. It will be displayed as pi_0 = exp(-rho)
    pprint(final_equation, use_unicode=False)

# Execute the function
find_pi_0()