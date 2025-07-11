import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    """
    k, c = sympy.symbols('k c')

    # Define jump sizes and their probabilities for large k
    jumps = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
        1:  sympy.Rational(1, 4) + c/k,
        2:  sympy.Rational(1, 4)
    }

    # Step 1: Calculate the expected drift mu_k
    # mu_k = E[X_{n+1} - k | X_n = k] = sum(jump * prob)
    mu_k = sum(jump * prob for jump, prob in jumps.items())
    mu_k = sympy.simplify(mu_k)
    
    print("Step 1: Calculating the one-step drift mu_k.")
    print(f"mu_k = (-2) * (1/4) + (-1) * (1/4 - c/k) + (1) * (1/4 + c/k) + (2) * (1/4)")
    print(f"mu_k = {mu_k}\n")

    # Step 2: Calculate the second moment (variance) of the jump sigma_k^2
    # sigma_k_sq = E[(X_{n+1} - k)^2 | X_n = k] = sum(jump^2 * prob)
    sigma_k_sq = sum(jump**2 * prob for jump, prob in jumps.items())
    sigma_k_sq = sympy.simplify(sigma_k_sq)
    
    print("Step 2: Calculating the one-step variance sigma_k^2.")
    print(f"sigma_k^2 = (-2)**2 * (1/4) + (-1)**2 * (1/4 - c/k) + (1)**2 * (1/4 + c/k) + (2)**2 * (1/4)")
    print(f"sigma_k^2 = {sigma_k_sq}\n")
    
    # Step 3: Apply the transience criterion
    # The chain is transient if 2 * lim_{k->inf} (k * mu_k) > lim_{k->inf} (sigma_k^2)
    # Let's calculate the limits
    limit_k_mu_k = sympy.limit(k * mu_k, k, sympy.oo)
    limit_sigma_k_sq = sympy.limit(sigma_k_sq, k, sympy.oo)
    
    print("Step 3: Applying the transience criterion.")
    print("The criterion for transience is: 2 * lim_{k->inf}(k * mu_k) > lim_{k->inf}(sigma_k^2)")
    print(f"The scaled asymptotic drift term, lim(k * mu_k), is: {limit_k_mu_k}")
    print(f"The asymptotic variance term, lim(sigma_k^2), is: {limit_sigma_k_sq}\n")
    
    # Step 4: Formulate and solve the inequality for c
    print("Step 4: Formulating the final inequality for c.")
    
    lhs = 2 * limit_k_mu_k
    rhs = limit_sigma_k_sq
    
    inequality = sympy.Gt(lhs, rhs)
    
    print("The inequality for transience is:")
    # We print each part of the equation separately as requested
    print(f"{lhs} > {rhs}")
    
    # Solve the inequality for c
    solution = sympy.solveset(inequality, c, domain=sympy.S.Reals)
    
    infimum = solution.inf
    
    print(f"\nSolving for c gives: {solution}")
    print(f"The set of values for c for which the chain is transient is {solution}.")
    print(f"The infimum of this set is {infimum}.")


if __name__ == '__main__':
    solve_markov_chain_transience()
