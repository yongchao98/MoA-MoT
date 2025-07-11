import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    The method involves analyzing the drift and variance of the chain's jumps.
    """
    
    # Define symbolic variables for c and k
    c, k = sympy.symbols('c k')

    # Step 1: Define the jump probabilities and calculate the expected drift (mu_k)
    print("Step 1: Calculate the expected drift mu_k = E[Jump]")
    
    # Jumps are -2, -1, +1, +2
    jumps = [-2, -1, 1, 2]
    
    # Probabilities for each jump from state k (for large k)
    probs = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
         1: sympy.Rational(1, 4) + c/k,
         2: sympy.Rational(1, 4)
    }

    mu_k = sum(jump * prob for jump, prob in probs.items())
    mu_k = sympy.simplify(mu_k)

    print(f"   mu_k = (-2)*({probs[-2]}) + (-1)*({probs[-1]}) + (1)*({probs[1]}) + (2)*({probs[2]})")
    print(f"   mu_k = {mu_k}\n")


    # Step 2: Calculate the limit of the variance of the jumps (sigma^2)
    print("Step 2: Calculate the asymptotic variance of the jumps sigma^2 = lim_{k->oo} E[Jump^2]")
    
    sigma_k_sq = sum(jump**2 * prob for jump, prob in probs.items())
    sigma_k_sq = sympy.simplify(sigma_k_sq)
    
    # As k -> infinity, terms with c/k vanish
    sigma_sq = sympy.limit(sigma_k_sq, k, sympy.oo)

    print(f"   sigma_k^2 = (-2)^2*({probs[-2]}) + (-1)^2*({probs[-1]}) + (1)^2*({probs[1]}) + (2)^2*({probs[2]})")
    print(f"   sigma_k^2 = {sigma_k_sq}")
    print(f"   Taking the limit as k -> oo, sigma^2 = {sigma_sq}\n")


    # Step 3: Apply the p-series test for transience
    print("Step 3: Apply the transience criterion")
    print("   The chain is transient if a corresponding p-series Sum(1/n^p) converges, which requires p > 1.")
    print("   The exponent p is determined by the drift and variance as: p = (2 * k * mu_k) / sigma^2")

    p = (2 * k * mu_k) / sigma_sq
    p = sympy.simplify(p)
    
    print(f"   p = (2 * k * ({mu_k})) / ({sigma_sq})")
    print(f"   p = {p}\n")

    # Step 4: Solve the inequality p > 1 for c
    print("Step 4: Solve the inequality for transience (p > 1)")
    
    # Using the final values to show the equation
    p_numerator = 8
    p_denominator = 5
    
    print(f"   The inequality is: ({p_numerator}*c) / {p_denominator} > 1")
    print(f"   Multiplying by {p_denominator}: {p_numerator}*c > {p_denominator}")
    print(f"   Dividing by {p_numerator}: c > {p_denominator}/{p_numerator}\n")


    # Step 5: Determine the infimum
    print("Step 5: Find the infimum")
    infimum_c = sympy.Rational(p_denominator, p_numerator)
    print(f"   The set of values for c that make the chain transient is ({infimum_c}, oo).")
    print(f"   The infimum (greatest lower bound) of this set is {infimum_c}.")

solve_markov_chain_transience()
<<<5/8>>>