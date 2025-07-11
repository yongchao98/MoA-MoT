import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    """
    # Define symbols for our variables
    c, k = sympy.symbols('c k')

    # Define the jumps (change in state) and their corresponding probabilities
    jumps_and_probs = {
        -2: sympy.Rational(1, 4),
        -1: sympy.Rational(1, 4) - c/k,
        1:  sympy.Rational(1, 4) + c/k,
        2:  sympy.Rational(1, 4)
    }

    # 1. Calculate the drift (mu_k)
    mu_k = sympy.simplify(sum(jump * prob for jump, prob in jumps_and_probs.items()))

    # 2. Calculate the second moment of the jump (M_2(k))
    m2_k = sympy.simplify(sum(jump**2 * prob for jump, prob in jumps_and_probs.items()))

    # For large k, the drift mu_k -> 0. The asymptotic variance sigma^2 is the limit of M_2(k).
    sigma_squared = sympy.limit(m2_k, k, sympy.oo)

    # 3. Apply the transience criterion.
    # The drift mu_k is of the form A/k. We find A by computing the limit of k * mu_k.
    A = sympy.limit(k * mu_k, k, sympy.oo)

    # The criterion for transience is A > sigma^2 / 2.
    # We will solve for the boundary case where A = sigma^2 / 2 to find the infimum.
    
    # We form the equation A = sigma^2 / 2
    boundary_equation = sympy.Eq(A, sigma_squared / 2)
    
    # 4. Solve for c to find the infimum value.
    infimum_c = sympy.solve(boundary_equation, c)

    # Print the derivation steps
    print(f"Step 1: Calculate the drift μ_k.")
    print(f"μ_k = E[X_n+1 - X_n | X_n = k] = {mu_k}")
    print("-" * 30)

    print(f"Step 2: Calculate the asymptotic variance of the jumps, σ².")
    print(f"The second moment of the jump is M₂(k) = {m2_k}.")
    print(f"The asymptotic variance σ² = lim(k->∞) M₂(k) = {sigma_squared}.")
    print("-" * 30)

    print("Step 3: Apply the transience criterion.")
    print("The drift μ_k is of the form A/k for large k.")
    print(f"A = lim(k->∞) k * μ_k = {A}.")
    print(f"The criterion for transience is A > σ²/2.")
    print("Substituting the derived values, the condition becomes:")
    
    # Print the final inequality with the numbers
    lhs = A
    rhs_numerator = sigma_squared.p
    rhs_denominator = sigma_squared.q * 2
    print(f"{lhs} > {rhs_numerator}/{rhs_denominator}")
    
    # Simplify and show the inequality for c
    inequality = sympy.Gt(A, sigma_squared / 2)
    print(f"Solving for c, we get: {inequality}")
    print("-" * 30)

    print("Step 4: Find the infimum.")
    print("The set of c for which the chain is transient is c > 5/8.")
    print(f"The infimum of this set is the boundary value, which is {infimum_c[0]}.")


solve_markov_chain_transience()
<<<5/8>>>