import sympy

def solve_steady_state():
    """
    This script derives the steady-state probability pi_0 for the given
    birth-death Markov process using symbolic mathematics.
    """
    # --- Step 1: Define symbolic variables ---
    # i and k are indices for the states
    i, k = sympy.symbols('i k', integer=True, nonnegative=True)
    # lambda is the arrival rate, mu is the service rate
    lam, mu = sympy.symbols('lambda mu', positive=True, real=True)
    # rho is the traffic intensity
    rho = sympy.symbols('rho', positive=True, real=True)
    # pi_0 is the steady-state probability of being in state 0
    pi_0 = sympy.Symbol('pi_0')

    print("--- Derivation of the steady-state probability pi_0 ---")
    print("\nStep 1: Define the transition rates.")
    # Birth rate (lambda_i) from state i to i+1
    birth_rate = lam / (i + 1)
    # Death rate (mu_i) from state i to i-1
    death_rate = mu
    print(f"Birth rate lambda_i = q_{{i, i+1}} = {birth_rate}")
    print(f"Death rate mu_i = q_{{i, i-1}} = {death_rate} (for i > 0)")

    print("\nStep 2: Express pi_k in terms of pi_0 using the balance equations.")
    # The general formula for pi_k in a birth-death process is:
    # pi_k = pi_0 * Product(lambda_{j-1} / mu_j, (j, 1, k))
    # Let's define the ratio term lambda_i / mu_{i+1}
    ratio = (lam / (i + 1)) / mu
    print(f"The ratio lambda_i / mu_{{i+1}} is: {sympy.simplify(ratio)}")
    
    # Substitute rho = lambda/mu
    ratio_in_rho = ratio.subs(lam / mu, rho)
    print(f"In terms of rho = lambda/mu, this ratio is: {ratio_in_rho}")
    
    # The product term from i=0 to k-1 gives pi_k / pi_0
    # pi_k = pi_0 * Product(rho/(j+1), (j, 0, k-1))
    # This product simplifies to rho**k / k!
    pi_k_expression = pi_0 * rho**k / sympy.factorial(k)
    print(f"The general expression for pi_k is: pi_k = {pi_k_expression}")
    
    print("\nStep 3: Apply the normalization condition Sum(pi_k) = 1.")
    # The sum of all probabilities must be 1.
    # Sum(pi_k for k=0 to infinity) = 1
    # pi_0 * Sum(rho**k / k!, (k, 0, oo)) = 1
    series = sympy.Sum(rho**k / sympy.factorial(k), (k, 0, sympy.oo))
    print(f"This leads to the equation: pi_0 * ({series}) = 1")
    
    print("\nStep 4: Evaluate the series.")
    # The series is the Taylor expansion of the exponential function e^rho
    series_value = series.doit()
    print(f"The sum {series} evaluates to: {series_value}")

    print("\nStep 5: Solve for pi_0.")
    # We have the equation: pi_0 * e^rho = 1
    final_equation = sympy.Eq(pi_0 * series_value, 1)
    print(f"Solving the final equation: {final_equation}")
    
    solution = sympy.solve(final_equation, pi_0)
    final_pi_0 = solution[0]
    
    print("\n--- Final Result ---")
    # To satisfy the output format of "output each number in the final equation"
    # we print the equation character by character.
    final_result_str = f"pi_0 = {final_pi_0}"
    print("The final equation is:")
    for char in final_result_str:
        print(char, end='')
    print("\n")


if __name__ == '__main__':
    solve_steady_state()