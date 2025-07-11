import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    """
    # Define symbolic variables for state k and parameter c
    k, c = sympy.symbols('k c')

    print("Step 1: Define the transition probabilities for large k.")
    # P(k, k+j) is the probability of moving from state k to k+j
    # Jumps are by -2, -1, +1, +2
    jumps = [-2, -1, 1, 2]
    probs = [
        sympy.Rational(1, 4),           # P(k, k-2)
        sympy.Rational(1, 4) - c/k,     # P(k, k-1)
        sympy.Rational(1, 4) + c/k,     # P(k, k+1)
        sympy.Rational(1, 4)            # P(k, k+2)
    ]
    
    print("Step 2: Calculate the drift (mu_k), which is the expected one-step change from state k.")
    # mu_k = E[X_{n+1} - X_n | X_n = k] = sum(jump * probability)
    mu_k = sum(jump * prob for jump, prob in zip(jumps, probs))
    mu_k = sympy.simplify(mu_k)
    print(f"The drift mu_k is: {mu_k}")

    print("\nStep 3: Calculate the second moment of the jump (sigma_k^2).")
    # sigma_k^2 = E[(X_{n+1} - k)^2 | X_n = k] = sum(jump^2 * probability)
    sigma_k_sq = sum(jump**2 * prob for jump, prob in zip(jumps, probs))
    sigma_k_sq = sympy.simplify(sigma_k_sq)
    print(f"The second moment sigma_k^2 is: {sigma_k_sq}")
    
    print("\nStep 4: Apply the transience criterion.")
    print("The critical boundary between transience and recurrence is given by the equation: mu_k = sigma_k^2 / (2*k)")
    
    # We want to solve mu_k = sigma_k^2 / (2*k) for c
    # (2*c/k) = (5/2) / (2*k)
    # 2*c = 5/4
    # This is an equation of the form A*c = B
    
    # Isolate the coefficient of c from the drift term
    A = sympy.simplify(mu_k * k / c)
    
    # Isolate the constant term from the right side of the criterion
    B = sympy.simplify((sigma_k_sq / (2*k)) * k)

    print("\nThe final equation to solve is of the form A*c = B, which simplifies from the criterion.")
    print(f"The number A (coefficient of c) is: {A}")
    print(f"The number B (constant term) is: {B}")

    # Solve for the critical value of c
    critical_c = B / A
    
    print(f"\nThe critical value for c is B/A = {B}/{A} = {critical_c}.")
    print("\nThe chain is transient for c > this critical value.")
    print("Therefore, the infimum of the set of c for which the chain is transient is this critical value.")

if __name__ == '__main__':
    solve_markov_chain_transience()
    # The final answer is 5/8
    # <<<5/8>>>