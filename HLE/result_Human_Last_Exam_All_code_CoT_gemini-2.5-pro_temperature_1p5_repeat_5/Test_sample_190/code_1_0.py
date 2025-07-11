import sympy

def solve_markov_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    """
    # Define symbols for state k and parameter c
    k = sympy.Symbol('k', positive=True, integer=True)
    c = sympy.Symbol('c')

    # Define the transition probabilities for large k as given
    P = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
        +1: sympy.Rational(1, 4) + c/k,
        +2: sympy.Rational(1, 4)
    }

    # Step 1: Calculate the local drift (mean increment) mu_k
    # mu_k = E[X_{n+1} - k | X_n = k] = sum(j * P_{k, k+j})
    mu_k = sum(jump * prob for jump, prob in P.items())
    mu_k = sympy.simplify(mu_k)

    print("Step 1: Calculate the local drift mu_k.")
    print(f"mu_k = (-2)*({P[-2]}) + (-1)*({P[-1]}) + (1)*({P[1]}) + (2)*({P[2]})")
    print(f"mu_k = {mu_k}\n")

    # The drift is of the form a/k. Let's find 'a'.
    a = (mu_k * k).expand()
    print(f"The drift is of the form a/k, where a = {a}.")
    print("-" * 40)

    # Step 2: Calculate the second moment of the increment, M2_k
    # M2_k = E[(X_{n+1} - k)^2 | X_n = k] = sum(j^2 * P_{k, k+j})
    M2_k = sum(jump**2 * prob for jump, prob in P.items())
    M2_k = sympy.simplify(M2_k)

    # Find its limit sigma_sq as k -> infinity
    sigma_sq = sympy.limit(M2_k, k, sympy.oo)
    
    print("Step 2: Calculate the limiting second moment of the increment, sigma^2.")
    print(f"M2_k = (-2)^2*({P[-2]}) + (-1)^2*({P[-1]}) + (1)^2*({P[1]}) + (2)^2*({P[2]})")
    print(f"M2_k = {M2_k}")
    print(f"The limit as k -> oo gives sigma^2 = {sigma_sq}.\n")
    print("-" * 40)
    
    # Step 3: Apply the criterion for transience: a > sigma^2 / 2
    print("Step 3: Apply the transience criterion a > sigma^2 / 2.")
    
    inequality_rhs = sigma_sq / 2
    
    print("The inequality for transience is:")
    print(f"{a} > {sigma_sq} / 2")
    print(f"{a} > {inequality_rhs}")

    # Step 4: Solve for c to find the set of transient values
    # We solve 2*c > 5/4
    infimum_c = sympy.solve(a > inequality_rhs, c)

    print(f"\nSolving for c, we find the condition for transience: {infimum_c}.")
    print("The set of c values for which the chain is transient is (5/8, infinity).")
    
    # The infimum is the lower bound of this interval.
    final_answer = inequality_rhs / 2 # from 2*c > 5/4 -> c > 5/8
    print(f"\nThe infimum of this set is {final_answer}.")

solve_markov_transience()