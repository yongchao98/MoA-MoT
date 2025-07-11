import sympy

def solve_steady_state_pi0():
    """
    This function symbolically derives the steady-state probability pi_0
    for the given birth-death process.
    """
    # Define the necessary symbolic variables
    pi_0 = sympy.Symbol('pi_0')
    rho = sympy.Symbol('rho', positive=True)  # rho = lambda/mu > 0
    k = sympy.Symbol('k', integer=True, nonnegative=True)

    # From the problem description, we derived that pi_k = pi_0 * (rho**k / k!)
    # The normalization condition is: Sum(pi_k for k=0 to inf) = 1
    # This can be written as: pi_0 * Sum(rho**k / k! for k=0 to inf) = 1

    # Define the term inside the summation
    sum_term = rho**k / sympy.factorial(k)

    # Create the infinite sum
    infinite_sum = sympy.Sum(sum_term, (k, 0, sympy.oo))

    # Evaluate the sum. Sympy recognizes this as the Taylor series for exp(rho).
    sum_value = infinite_sum.doit()

    # Set up the normalization equation: pi_0 * sum_value = 1
    normalization_eq = sympy.Eq(pi_0 * sum_value, 1)

    # Solve the equation for pi_0
    pi_0_solution = sympy.solve(normalization_eq, pi_0)

    # The solution is a list, we take the first element
    final_expression = pi_0_solution[0]
    
    # Create the final equation object for printing
    final_equation = sympy.Eq(pi_0, final_expression)
    
    print("The relationship between pi_k and pi_0 is: pi_k = pi_0 * rho**k / k!")
    print("The normalization condition is: pi_0 * Sum(rho**k / k!) = 1")
    print(f"The infinite sum evaluates to: {sum_value}")
    print(f"The equation becomes: {pi_0} * {sum_value} = 1")
    print("\nSolving for pi_0, we get the final equation:")
    
    # To satisfy the prompt "output each number in the final equation",
    # we will print the equation in a structured way. The expression is exp(-1 * rho).
    print(f"{final_equation.lhs} = {sympy.latex(final_equation.rhs)}")
    
solve_steady_state_pi0()