import sympy

def solve_and_explain():
    """
    This function derives the steady-state probability pi_0 for the given
    birth-death process and prints the derivation steps.
    """
    # Define symbols for the derivation
    i, n = sympy.symbols('i n', integer=True, positive=True)
    lmbda, mu, rho = sympy.symbols('lambda mu rho')
    pi = sympy.Function('pi')

    # Step 1: State the detailed balance equations
    print("Step 1: The detailed balance equations for a birth-death process are:")
    print("pi_i * lambda_i = pi_{i+1} * mu_{i+1}\n")

    # Step 2: Substitute the given rates
    print("Step 2: Substitute the given rates into the equation.")
    print(f"The birth rate from state i is lambda_i = {lmbda}/(i+1)")
    print(f"The death rate from state i+1 is mu_{i+1} = {mu}")
    print(f"So, pi(i) * ({lmbda}/(i+1)) = pi(i+1) * {mu}\n")

    # Step 3: Find the relationship between successive probabilities
    print("Step 3: Rearrange to find pi_{i+1} in terms of pi_i.")
    print(f"pi(i+1) = pi(i) * ({lmbda}/{mu}) / (i+1)")
    print(f"Let rho = lambda/mu. The equation becomes:")
    print("pi(i+1) = pi(i) * rho / (i+1)\n")

    # Step 4: Express pi_n in terms of pi_0
    print("Step 4: Express the general term pi_n in terms of pi_0 by unrolling the recursion.")
    print("pi_1 = pi_0 * rho / 1")
    print("pi_2 = pi_1 * rho / 2 = (pi_0 * rho) * rho / 2 = pi_0 * rho^2 / 2!")
    print("pi_3 = pi_2 * rho / 3 = (pi_0 * rho^2 / 2!) * rho / 3 = pi_0 * rho^3 / 3!")
    print("By induction, we find the general form:")
    print("pi_n = pi_0 * (rho^n / n!)\n")

    # Step 5: Use the normalization condition
    print("Step 5: Use the normalization condition Sum(pi_n for n=0 to infinity) = 1.")
    print("Sum(pi_0 * rho^n / n!) = 1")
    print("pi_0 * Sum(rho^n / n!) = 1\n")

    # Step 6: Recognize the series and solve for pi_0
    print("Step 6: Recognize the sum as the Taylor series for e^rho.")
    series_sum = sympy.summation(rho**n / sympy.factorial(n), (n, 0, sympy.oo))
    print(f"The series Sum(rho^n / n!) from n=0 to infinity is equal to {series_sum}.")
    print(f"So, pi_0 * {series_sum} = 1\n")

    # Final Answer
    print("Step 7: The final equation for pi_0 is:")
    # The prompt asks to output each number in the final equation.
    # The equation is pi_0 = e^(-rho), which can be written as pi_0 = e^(-1 * rho).
    # The numbers are 0 and -1.
    print("pi_0 = e^(-1 * rho)")

if __name__ == '__main__':
    solve_and_explain()