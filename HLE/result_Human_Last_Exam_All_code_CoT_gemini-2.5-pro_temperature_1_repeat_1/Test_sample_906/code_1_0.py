import sympy as sp

def solve_steady_state_pi0():
    """
    Solves for the steady-state probability pi_0 for the given birth-death process.
    The relationship between pi_n and pi_0 is pi_n = pi_0 * rho^n / n!
    The sum of all probabilities must be 1: Sum(pi_n for n=0 to inf) = 1.
    This means pi_0 * Sum(rho^n / n! for n=0 to inf) = 1.
    We solve for pi_0.
    """
    # Define the symbolic variables
    n = sp.symbols('n', integer=True, nonnegative=True)
    rho = sp.symbols('rho', real=True, positive=True)
    pi_0 = sp.symbols('pi_0')

    # The term inside the summation for pi_n / pi_0
    term = rho**n / sp.factorial(n)

    # Calculate the infinite sum using sympy's summation function
    # This sum represents Sum(pi_n / pi_0)
    total_sum = sp.summation(term, (n, 0, sp.oo))

    # The normalization equation is pi_0 * total_sum = 1
    normalization_eq = sp.Eq(pi_0 * total_sum, 1)
    
    # Solve the equation for pi_0
    solution = sp.solve(normalization_eq, pi_0)
    
    # The solution is a list, get the first element
    pi_0_expression = solution[0]
    
    # Print the final equation
    print(f"The normalization equation is: {pi_0} * ({total_sum}) = 1")
    print(f"The expression for the steady-state probability pi_0 is:")
    print(f"pi_0 = {pi_0_expression}")

if __name__ == "__main__":
    solve_steady_state_pi0()