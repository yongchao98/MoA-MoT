import sympy

def find_transience_infimum():
    """
    Calculates the infimum of c for which the given Markov chain is transient.

    This script follows a standard method for analyzing one-dimensional random walks:
    1. Define the jump probabilities symbolically.
    2. Compute the asymptotic mean drift (mu_k) and jump variance (sigma_k^2).
    3. Use a diffusion approximation (Lamperti's criterion) to establish the condition for transience.
    4. Solve the resulting inequality for 'c' to find the desired infimum.
    """
    # Step 1: Define symbols and transition probabilities for large k.
    print("Step 1: Defining the symbolic representation of the problem.")
    k, c = sympy.symbols('k c')

    # Transition probabilities P_{k, k+j} for a jump of size j
    probabilities = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
         1: sympy.Rational(1, 4) + c/k,
         2: sympy.Rational(1, 4)
    }
    print("Transition probabilities defined.")
    print("-" * 40)

    # Step 2: Calculate the mean drift mu_k.
    print("Step 2: Calculating the mean drift mu_k = E[X_{n+1} - k | X_n = k].")
    mu_k = sum(jump * prob for jump, prob in probabilities.items())
    mu_k = sympy.simplify(mu_k)
    print(f"Calculated mean drift mu_k = {mu_k}")
    print("-" * 40)

    # Step 3: Calculate the second moment and asymptotic variance.
    print("Step 3: Calculating the second moment B_k = E[(X_{n+1} - k)^2 | X_n = k].")
    B_k = sum(jump**2 * prob for jump, prob in probabilities.items())
    B_k = sympy.simplify(B_k)
    print(f"Calculated second moment B_k = {B_k}")

    print("\nFor large k, the variance sigma_k^2 = B_k - mu_k^2 approaches B_k because mu_k -> 0.")
    sigma_sq_limit = B_k
    print(f"The asymptotic variance sigma^2 is {sigma_sq_limit}.")
    print("-" * 40)

    # Step 4: Apply the criterion for transience.
    print("Step 4: Applying Lamperti's criterion for transience.")
    print("The chain is transient if the exponent 'p' from the analysis is greater than 1.")
    
    # The exponent 'p' is the coefficient of 1/k in the term 2*mu_k / sigma_k^2.
    criterion_term = 2 * mu_k / sigma_sq_limit
    p = criterion_term.coeff(1/k)
    
    print(f"The key term 2*mu_k / sigma_k^2 simplifies to: {criterion_term}")
    print(f"The resulting exponent 'p' is: {p}")
    
    p_num, p_den = p.as_numer_denom()
    p_c_coeff = p_num.coeff(c)

    # Step 5: Formulate and solve the inequality for c.
    print("\nThe condition for transience is p > 1.")
    print("This gives the following inequality for c:")
    # We output the numbers in the equation as requested
    print(f"Equation: {p_c_coeff} * c / {p_den} > 1")

    # Solving for c gives c > p_den / p_c_coeff
    solution_threshold = p_den / p_c_coeff
    print(f"Solving for c, we get: c > {p_den}/{p_c_coeff}")
    print("-" * 40)

    # Step 6: State the final answer.
    print("Step 6: Determining the infimum.")
    print(f"The Markov chain is transient for c in the set ({solution_threshold}, infinity).")
    print(f"The infimum of this set is the boundary value.")
    print(f"\nFinal Answer: The infimum is {solution_threshold}.")

if __name__ == "__main__":
    find_transience_infimum()