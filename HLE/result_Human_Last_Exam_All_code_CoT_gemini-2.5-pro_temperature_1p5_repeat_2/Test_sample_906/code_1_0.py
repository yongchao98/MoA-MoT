import sympy

def solve_steady_state_pi0():
    """
    This function symbolically derives the steady-state probability pi_0
    for the given birth-death process.
    """
    # Define the symbols we will use
    pi_0 = sympy.Symbol('pi_0')
    rho = sympy.Symbol('rho', positive=True)
    n = sympy.Symbol('n', integer=True, nonnegative=True)

    # As derived in the explanation, pi_n = pi_0 * (rho**n / n!)
    # The sum of all probabilities must be 1: Sum(pi_n for n=0 to inf) = 1
    # This simplifies to: pi_0 * Sum(rho**n / n!) = 1

    # The sum is the Taylor series for exp(rho)
    series_sum = sympy.Sum(rho**n / sympy.factorial(n), (n, 0, sympy.oo)).doit()

    # We now have the normalization equation: pi_0 * exp(rho) = 1
    # Let's represent the terms of this final equation
    # LHS: pi_0 * series_sum
    # RHS: 1
    lhs = pi_0 * series_sum
    rhs = 1
    
    final_equation = sympy.Eq(lhs, rhs)
    
    print("The final normalization equation is:")
    # To make it more readable, let's substitute the sum with its simplified form in the printout
    print(f"{pi_0} * exp({rho}) = 1")
    
    # Solve the equation for pi_0
    solution = sympy.solve(final_equation, pi_0)

    print("\nThe solution for pi_0 is:")
    # The result of solve is a list, so we take the first element
    print(solution[0])

if __name__ == '__main__':
    solve_steady_state_pi0()