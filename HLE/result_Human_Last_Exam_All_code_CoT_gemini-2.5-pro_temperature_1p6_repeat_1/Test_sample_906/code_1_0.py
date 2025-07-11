import sympy

def derive_pi0():
    """
    Derives and prints the steady-state probability pi_0 for the given birth-death process.
    The derivation is presented step-by-step.
    """

    # --- Symbolic setup for the derivation ---
    pi_0_sym = sympy.Symbol('pi_0')
    rho_sym = sympy.Symbol('rho')
    k_sym = sympy.Symbol('k', integer=True, nonnegative=True)

    # --- Derivation Steps ---
    print("Derivation of the steady-state probability pi_0:")
    print("==================================================")

    # Step 1: Express pi_k in terms of pi_0
    print("Step 1: Express pi_k in terms of pi_0.")
    print("For a general birth-death process, pi_k = pi_0 * product_{i=0 to k-1} (lambda_i / mu_{i+1}).")
    print("Given the rates lambda_i = lambda / (i+1) and mu_i = mu, the ratio becomes:")
    print("lambda_i / mu_{i+1} = (lambda / (i+1)) / mu = (lambda/mu) / (i+1)")
    print("\nThe product term is product_{i=0 to k-1} [ (lambda/mu) / (i+1) ] = (lambda/mu)^k * (1/k!)")
    print(f"Letting rho = lambda/mu, we get the simplified relation: pi_k = pi_0 * (rho^k / k!)\n")

    # Step 2: Apply the normalization condition
    print("Step 2: Use the normalization condition Sum(pi_k for k=0 to infinity) = 1.")
    print("Sum_{k=0 to infinity} [pi_0 * (rho^k / k!)] = 1")
    print("This can be rewritten as: pi_0 * Sum_{k=0 to infinity} [rho^k / k!] = 1\n")

    # Step 3: Identify the series and solve
    print("Step 3: Recognize the sum as the Taylor series for e^rho and solve.")
    sum_expr = sympy.Sum(rho_sym**k_sym / sympy.factorial(k_sym), (k_sym, 0, sympy.oo))
    sum_val = sum_expr.doit()
    print(f"The sum Sum_{k=0 to infinity} [rho^k / k!] is the series for {sum_val}.")
    print(f"So, the equation becomes: pi_0 * {sum_val} = 1\n")

    # Step 4: State the final result
    print("Step 4: The final equation for pi_0.")
    print("Solving for pi_0 gives: pi_0 = 1 / e^rho = e^(-rho)")

    # Output the components of the final equation as requested
    print("\nFinal Equation Breakdown:")
    print("pi_0 = e**(-1 * rho)")
    print("-------------------------")
    print("Symbol      | Value/Meaning")
    print("-------------------------")
    print("e           | Euler's number (base of natural logarithm)")
    print("-1          | Coefficient in the exponent")
    print("rho         | The traffic intensity (lambda / mu)")
    print("-------------------------")

if __name__ == "__main__":
    derive_pi0()