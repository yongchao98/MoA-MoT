import sympy

def solve_steady_state():
    """
    Symbolically derives the steady-state probability pi_0 for the given
    birth-death process.
    """
    # Define symbolic variables for the derivation
    pi_0 = sympy.Symbol('pi_0')
    rho = sympy.Symbol('rho', positive=True) # rho = lambda/mu, must be positive
    n = sympy.Symbol('n', integer=True, nonnegative=True)

    # From the detailed balance equations, we derive pi_n = pi_0 * rho**n / n!
    pi_n = pi_0 * rho**n / sympy.factorial(n)

    # The sum of all probabilities must equal 1 (Normalization Condition)
    # We define the sum from n=0 to infinity
    # Sum(pi_n) for n=0..oo
    normalization_sum = sympy.Sum(pi_n, (n, 0, sympy.oo))

    # Evaluate the sum. sympy knows the Taylor series for exp(rho)
    evaluated_sum = normalization_sum.doit()

    # The normalization equation is Sum(pi_n) = 1
    # which becomes: pi_0 * exp(rho) = 1
    normalization_equation = sympy.Eq(evaluated_sum, 1)

    print("The derived relationship between probabilities is: pi_n = pi_0 * (rho^n / n!)")
    print("The normalization condition requires the sum of all probabilities to be 1:")
    print(f"Sum({pi_n}) = 1")
    print("\nThis leads to the equation:")
    # The final equation has one number: 1
    print(f"{normalization_equation.lhs} = {normalization_equation.rhs}")


    # Solve the equation for pi_0
    solution = sympy.solve(normalization_equation, pi_0)

    # The result is a list, so we extract the first element
    pi_0_solution = solution[0]

    print("\nSolving for pi_0 gives the final expression:")
    print(f"pi_0 = {pi_0_solution}")

solve_steady_state()