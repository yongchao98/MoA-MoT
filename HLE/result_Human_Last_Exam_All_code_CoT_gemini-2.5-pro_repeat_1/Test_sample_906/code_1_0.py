import sympy

def solve_steady_state():
    """
    Derives the steady-state probability pi_0 for the given birth-death process.
    """
    # Define symbols for the derivation
    k = sympy.Symbol('k', integer=True, nonnegative=True)
    rho = sympy.Symbol('rho')
    pi_0 = sympy.Symbol('pi_0')

    print("Step 1: The general formula for steady-state probabilities in a birth-death process is:")
    print("pi_k = pi_0 * product_{i=0 to k-1} (lambda_i / mu_{i+1})\n")

    print("Step 2: We substitute the given rates:")
    print("Birth rate lambda_i = lambda / (i+1)")
    print("Death rate mu_i = mu")
    print("The ratio is lambda_i / mu_{i+1} = (lambda / (i+1)) / mu = (lambda/mu) / (i+1).")
    print("Let rho = lambda/mu. The ratio becomes rho / (i+1).\n")

    print("Step 3: We compute the product term to find the relationship between pi_k and pi_0.")
    # Product is (rho/1) * (rho/2) * ... * (rho/k) = rho^k / k!
    pi_k_expression = rho**k / sympy.factorial(k)
    print("product_{i=0 to k-1} (rho / (i+1)) = (rho/1)*(rho/2)*...*(rho/k) = rho**k / k!")
    print(f"So, pi_k = pi_0 * {pi_k_expression}\n")

    print("Step 4: We use the normalization condition that the sum of all probabilities is 1:")
    print("sum_{k=0 to infinity} pi_k = 1")
    print("This leads to: pi_0 * (sum_{k=0 to infinity} (rho**k / k!)) = 1\n")

    print("Step 5: We recognize the infinite sum.")
    # The sum is the Taylor series for e^rho
    series_sum = sympy.summation(rho**k / sympy.factorial(k), (k, 0, sympy.oo))
    print(f"The series sum_{k=0 to infinity} (rho**k / k!) is the Taylor expansion for e**rho.")
    print(f"In sympy, this evaluates to: {series_sum}\n")

    print("Step 6: We solve for pi_0.")
    print(f"Substituting the sum back, we get: pi_0 * {series_sum} = 1")
    final_pi_0_expr = 1 / series_sum
    print(f"Solving for pi_0 gives: pi_0 = 1 / {series_sum}\n")

    print("--- FINAL EQUATION ---")
    # To satisfy "output each number in the final equation", we'll be explicit.
    # The equation is pi_0 = e^(-rho), which can be seen as e**(-1 * rho)
    base = 'e'
    exponent_coeff = -1
    variable = 'rho'
    print(f"The final expression for the steady-state probability pi_0 is:")
    print(f"pi_0 = {base}**({exponent_coeff} * {variable})")
    print(f"pi_0 = {final_pi_0_expr}")

solve_steady_state()