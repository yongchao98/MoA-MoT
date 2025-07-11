import sympy
from sympy import symbols, summation, oo, factorial, Eq, solve, exp

def solve_steady_state():
    """
    This function uses sympy to derive the steady-state probability pi_0
    for the given birth-death process.
    """
    # Define the symbolic variables
    # rho represents lambda/mu
    # k is the index for the state and the summation
    # pi_0 is the steady-state probability of being in state 0
    rho = symbols('rho', positive=True)
    k = symbols('k', integer=True, nonnegative=True)
    pi_0 = symbols('pi_0')

    # The relationship between pi_k and pi_0 is pi_k = pi_0 * (rho**k / k!)
    # The normalization condition is Sum(pi_k for k>=0) = 1, which means
    # pi_0 * Sum(rho**k / k! for k>=0) = 1.

    # We first need to evaluate the sum of the series.
    series_term = rho**k / factorial(k)
    series_sum = summation(series_term, (k, 0, oo))

    # Now, we form the normalization equation.
    # The left side is pi_0 * the sum we just calculated, and the right side is 1.
    equation = Eq(pi_0 * series_sum, 1)

    # Solve the equation for pi_0. The result is a list with one element.
    solution = solve(equation, pi_0)[0]
    
    # Print the results in a step-by-step manner.
    print("The relationship between steady-state probabilities is: pi_k = pi_0 * (rho^k / k!)")
    print("The normalization condition is: pi_0 * Sum(rho^k / k! for k=0 to infinity) = 1")
    print(f"\nThe sum of the series Sum(rho^k / k!) is: {series_sum}")
    
    # Display the final equation using the calculated sum
    final_equation_str = f"pi_0 * ({series_sum}) = 1"
    print(f"\nThis gives the final equation: {final_equation_str}")
    
    # Display the solution
    print(f"\nSolving for pi_0, we find:")
    print(f"pi_0 = {solution}")

# Execute the function to find and print the solution
solve_steady_state()