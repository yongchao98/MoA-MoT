import sympy

def solve_steady_state():
    """
    This function derives the steady-state probability pi_0 for the given
    birth-death process using symbolic mathematics.
    """
    # Define the symbols we will use
    rho = sympy.Symbol('rho')
    k = sympy.Symbol('k', integer=True, nonnegative=True)
    pi_0 = sympy.Symbol('pi_0')

    # The relationship between pi_k and pi_0 is pi_k = pi_0 * rho**k / k!
    # The normalization condition is Sum(pi_k for k=0 to infinity) = 1.
    # This leads to the equation: pi_0 * Sum(rho**k / k!) = 1.

    # Use sympy to compute the infinite sum part of the equation.
    # The sum is Sum(rho**k / k! for k from 0 to infinity)
    try:
        infinite_sum = sympy.summation(rho**k / sympy.factorial(k), (k, 0, sympy.oo))
    except Exception as e:
        print(f"Could not compute the sum symbolically: {e}")
        return

    # Now we have the equation: pi_0 * infinite_sum = 1
    # We can represent this equation in sympy
    equation = sympy.Eq(pi_0 * infinite_sum, 1)

    # Solve the equation for pi_0
    solution = sympy.solve(equation, pi_0)

    # The solution is a list, so we extract the first element.
    final_expression = solution[0]

    # Print the derivation steps and the final result
    print("Step 1: The relationship between probabilities is pi_k = pi_0 * (rho**k / k!)")
    print("\nStep 2: The normalization condition is Sum(pi_k) = 1, which means:")
    print(f"pi_0 * Sum(rho**k / k!) = 1")
    print(f"\nStep 3: The sum Sum(rho**k / k!) is the Taylor series for e**rho, which is: {infinite_sum}")
    print(f"\nStep 4: The equation becomes pi_0 * {infinite_sum} = 1")
    print("\nFinal Answer: Solving for pi_0, we get:")
    print(f"pi_0 = {final_expression}")

solve_steady_state()