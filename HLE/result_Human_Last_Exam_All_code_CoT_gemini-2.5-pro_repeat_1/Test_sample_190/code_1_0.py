import sympy

def solve_markov_transience():
    """
    This function calculates the infimum of c for which the given Markov chain is transient.
    """
    # Define symbolic variables for state k and parameter c
    k, c = sympy.symbols('k c')

    # Step 1: Define jump sizes and their probabilities for large k
    jumps = [-2, 2, -1, 1]
    probabilities = [
        sympy.Rational(1, 4), 
        sympy.Rational(1, 4), 
        sympy.Rational(1, 4) - c/k, 
        sympy.Rational(1, 4) + c/k
    ]

    # Step 2: Calculate the mean drift mu_k
    mu_k = sum(jump * prob for jump, prob in zip(jumps, probabilities))
    mu_k = sympy.simplify(mu_k)
    
    print(f"The mean drift at state k is: mu_k = {mu_k}")

    # Determine the asymptotic coefficient A, where mu_k ~ A/k
    A = sympy.limit(mu_k * k, k, sympy.oo)
    print(f"The drift is asymptotically A/k, with A = {A}")
    print("-" * 30)

    # Step 3: Calculate the second moment M2_k and asymptotic variance sigma^2
    M2_k = sum(jump**2 * prob for jump, prob in zip(jumps, probabilities))
    M2_k = sympy.simplify(M2_k)

    # For large k, the variance sigma_k^2 approaches M2_k since mu_k -> 0.
    sigma_sq = sympy.limit(M2_k, k, sympy.oo)
    
    print(f"The second moment of the jump is: M2_k = {M2_k}")
    print(f"The asymptotic variance of the jump is: sigma^2 = {sigma_sq}")
    print("-" * 30)

    # Step 4: Apply the transience criterion 2*A / sigma^2 > 1 and solve for c
    print("The criterion for the chain to be transient is: 2*A / sigma^2 > 1")
    
    # Create the inequality
    inequality_lhs = 2 * A / sigma_sq
    print(f"Substituting A and sigma^2: {inequality_lhs} > 1")

    # The inequality simplifies to 8*c / 5 > 1
    # We output the numbers in the final simplified form of the inequality
    num_8 = 8
    num_5 = 5
    print(f"Solving for c, the inequality becomes: {num_8}*c > {num_5}")
    
    # Step 5: Determine the infimum
    infimum = sympy.Rational(num_5, num_8)
    print(f"\nThis implies that the chain is transient for c > {infimum}.")
    print(f"The set of c values for transience is ({infimum}, infinity).")
    print(f"Therefore, the infimum of this set is {infimum}.")

solve_markov_transience()
<<<5/8>>>