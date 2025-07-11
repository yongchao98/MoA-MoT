import sympy

def solve_markov_chain_transience():
    """
    This function calculates the infimum of c for which the given Markov chain is transient.
    It follows a standard procedure based on the drift and variance of the chain.
    """
    
    # Define symbolic variables for k and c
    k = sympy.Symbol('k')
    c = sympy.Symbol('c')

    # Define the transition probabilities for large k as given in the problem
    P = {
        -2: sympy.Rational(1, 4),           # P_{k, k-2}
        -1: sympy.Rational(1, 4) - c/k,     # P_{k, k-1}
        1:  sympy.Rational(1, 4) + c/k,     # P_{k, k+1}
        2:  sympy.Rational(1, 4)            # P_{k, k+2}
    }
    
    # --- Step 1: Calculate the mean drift mu_k ---
    # The mean drift is the expected jump size from state k.
    # mu_k = E[X_{n+1} - X_n | X_n = k]
    mu_k = sum(jump * prob for jump, prob in P.items())
    mu_k = sympy.simplify(mu_k)
    
    print("Step 1: Calculate the mean drift mu_k.")
    print(f"mu_k = (-2) * ({P[-2]}) + (-1) * ({P[-1]}) + (1) * ({P[1]}) + (2) * ({P[2]})")
    print(f"After simplification, mu_k = {mu_k}\n")
    
    # --- Step 2: Identify the parameter alpha ---
    # According to the transience criterion, the drift mu_k should be of the form alpha/k.
    # We find alpha by calculating the limit of k * mu_k as k -> infinity.
    alpha = sympy.limit(k * mu_k, k, sympy.oo)
    
    print("Step 2: Identify the parameter alpha from the drift equation mu_k = alpha/k.")
    print(f"alpha = lim_{{k->oo}} (k * mu_k) = lim_{{k->oo}} (k * ({mu_k}))")
    print(f"alpha = {alpha}\n")

    # --- Step 3: Calculate the asymptotic variance sigma^2 ---
    # First, we calculate the second moment of the jump size from state k.
    second_moment_k = sum(jump**2 * prob for jump, prob in P.items())
    second_moment_k = sympy.simplify(second_moment_k)
    
    # The variance is Var = E[X^2] - (E[X])^2.
    # As k -> oo, the drift mu_k -> 0.
    # So, the limit of the variance is the limit of the second moment.
    sigma_squared = sympy.limit(second_moment_k, k, sympy.oo)

    print("Step 3: Calculate the asymptotic variance sigma^2.")
    print(f"Second Moment E[(Jump)^2] = (-2)^2*({P[-2]}) + (-1)^2*({P[-1]}) + (1)^2*({P[1]}) + (2)^2*({P[2]})")
    print(f"Second Moment = {second_moment_k}")
    print(f"sigma^2 = lim_{{k->oo}} (Second Moment) = {sigma_squared}\n")

    # --- Step 4: Apply the transience criterion ---
    # The Markov chain is transient if alpha > sigma^2 / 2.
    # We construct and solve this inequality for c.
    print("Step 4: Apply the transience criterion alpha > sigma^2 / 2 and solve for c.")
    
    # Get the coefficient of c from alpha, which is a linear expression in c.
    # For alpha = 2*c, the coefficient is 2.
    alpha_coeff = sympy.Poly(alpha, c).coeffs()[0]
    
    # Calculate the value for the right side of the inequality.
    sigma_squared_div_2 = sigma_squared / 2

    # Print the final equation with all numbers substituted
    print("The transience condition is:")
    print(f"     alpha > sigma^2 / 2")
    print("Substituting the calculated values:")
    print(f"     {alpha} > {sigma_squared} / 2")
    print(f"     {alpha} > {sigma_squared_div_2}")
    print(f"Solving for c:")
    print(f"     c > {sigma_squared_div_2} / {alpha_coeff}")
    
    infimum_c = sigma_squared_div_2 / alpha_coeff
    print(f"     c > {infimum_c}\n")

    print(f"The set of values of c for which the chain is transient is ({infimum_c}, infinity).")
    print(f"The infimum of this set is {infimum_c}.")

if __name__ == '__main__':
    solve_markov_chain_transience()