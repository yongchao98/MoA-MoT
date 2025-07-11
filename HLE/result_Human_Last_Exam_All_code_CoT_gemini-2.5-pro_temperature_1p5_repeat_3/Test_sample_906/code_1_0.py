import sympy as sp

def solve_steady_state():
    """
    This function derives the steady-state probability pi_0 for the given
    birth-death process using symbolic mathematics.
    """
    # 1. Define symbolic variables
    pi_0 = sp.Symbol('pi_0')
    rho = sp.Symbol('rho')
    i = sp.Symbol('i', integer=True, nonnegative=True)

    print("Step 1: Define the relationship between probabilities.")
    print("For a birth-death process, the detailed balance equation is pi_i * lambda_i = pi_{i+1} * mu_{i+1}.")
    print("Given the rates lambda_i = lambda / (i+1) and mu_i = mu, we have:")
    print("pi_i * (lambda / (i+1)) = pi_{i+1} * mu")
    print("Rearranging for pi_{i+1} and letting rho = lambda / mu gives:")
    print("pi_{i+1} = pi_i * (rho / (i+1))")
    print("\nBy applying this rule repeatedly, we find the general expression for pi_i in terms of pi_0:")
    # General expression for pi_i
    pi_i_expr = pi_0 * rho**i / sp.factorial(i)
    print(f"pi_i = pi_0 * rho**i / i!  which is  {pi_i_expr}")
    print("-" * 50)

    print("Step 2: Apply the normalization condition.")
    print("The sum of all probabilities must equal 1:")
    print("Sum(pi_i for i from 0 to infinity) = 1")
    sum_expr = sp.Sum(pi_i_expr, (i, 0, sp.oo))
    print("Substituting the expression for pi_i, the equation is:")
    print(f"{sum_expr} = 1")
    print("-" * 50)

    print("Step 3: Solve the equation for pi_0.")
    print("We can factor pi_0 out of the summation:")
    factored_sum = pi_0 * sp.Sum(rho**i / sp.factorial(i), (i, 0, sp.oo))
    print(f"{factored_sum} = 1")
    print("\nThe summation is the well-known Taylor series expansion for the exponential function e**rho.")
    
    # Evaluate the sum
    sum_value = factored_sum.doit()
    print(f"So, the equation becomes:")
    final_eq = sp.Eq(sum_value, 1)
    print(f"{sp.pretty(final_eq, use_unicode=False)}")

    # Solve for pi_0
    solution = sp.solve(final_eq, pi_0)
    
    print("\nSolving for pi_0, we get:")
    final_solution_eq = sp.Eq(pi_0, solution[0])
    print(f"{sp.pretty(final_solution_eq, use_unicode=False)}")
    print("-" * 50)

if __name__ == '__main__':
    solve_steady_state()