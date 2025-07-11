import sympy as sp

def solve_steady_state():
    """
    Derives the steady-state probability pi_0 for the given birth-death process
    using symbolic mathematics.
    """
    # Define the symbolic variables
    pi_0 = sp.Symbol('pi_0')
    rho = sp.Symbol('rho')
    k = sp.Symbol('k', integer=True, nonnegative=True)

    # The expression for pi_k in terms of pi_0 is pi_k = pi_0 * rho**k / k!
    pi_k_expr = pi_0 * rho**k / sp.factorial(k)

    # The normalization condition is Sum(pi_k) from k=0 to infinity = 1
    # We can write this as pi_0 * Sum(rho**k / k!) = 1
    sum_expr = sp.Sum(rho**k / sp.factorial(k), (k, 0, sp.oo))
    
    # Evaluate the sum. Sympy recognizes this as the Taylor series for exp(rho).
    evaluated_sum = sum_expr.doit()

    # Form the final equation before solving for pi_0
    final_equation = sp.Eq(pi_0 * evaluated_sum, 1)

    # Solve the equation for pi_0
    solution = sp.solve(final_equation, pi_0)
    
    # The result is the first (and only) element in the solution list
    pi_0_solution = solution[0]

    # Print the derivation steps
    print("The steady-state probability pi_k is given by: pi_0 * rho**k / k!")
    print("The normalization condition requires the sum of all probabilities to be 1:")
    print("pi_0 * Sum(rho**k / k!, (k, 0, oo)) = 1\n")
    print("The sum is the Taylor series for exp(rho), so the final equation is:")
    
    # Print the terms of the final equation: pi_0 * exp(rho) = 1
    lhs_coeff = final_equation.lhs.coeff(pi_0)
    rhs_const = final_equation.rhs
    
    # sp.pretty_print is used for a nicer mathematical layout
    print(f"pi_0 * ({sp.pretty(lhs_coeff)}) = {sp.pretty(rhs_const)}\n")
    
    print("Solving for pi_0, we get:")
    print(f"pi_0 = {sp.pretty(pi_0_solution)}")


if __name__ == "__main__":
    solve_steady_state()