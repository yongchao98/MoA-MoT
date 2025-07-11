import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for a transient Markov chain
    using the Lamperti-Kesten criterion.
    """
    # 1. Define symbolic variables
    c, k = sympy.symbols('c k')

    # 2. Define transition probabilities for large k
    p_km2 = sympy.Rational(1, 4)  # P(k, k-2)
    p_kp2 = sympy.Rational(1, 4)  # P(k, k+2)
    p_km1 = sympy.Rational(1, 4) - c/k  # P(k, k-1)
    p_kp1 = sympy.Rational(1, 4) + c/k  # P(k, k+1)

    print("Step 1: Calculate the expected drift mu_k = E[X_{n+1} - X_n | X_n = k].")
    
    # Calculate mu_k, the expected drift from state k
    jumps = [-2, -1, 1, 2]
    probs = [p_km2, p_km1, p_kp1, p_kp2]
    
    mu_k = sympy.simplify(sum(j * p for j, p in zip(jumps, probs)))
    
    print(f"The expected drift is mu_k = (-2) * ({p_km2}) + (-1) * ({p_km1}) + (1) * ({p_kp1}) + (2) * ({p_kp2})")
    print(f"mu_k = {mu_k}\n")
    
    print("Step 2: Calculate the asymptotic variance of the jump, sigma^2.")
    
    # Calculate E[Jump^2] and then the limit for the variance
    E_sq_jump_k = sympy.simplify(sum((j**2) * p for j, p in zip(jumps, probs)))
    # Since mu_k -> 0 as k -> oo, the asymptotic variance sigma^2 is the limit of E[Jump^2]
    sigma_sq = sympy.limit(E_sq_jump_k, k, sympy.oo)

    print(f"The expected squared jump is E[Jump^2] = (-2)^2 * ({p_km2}) + (-1)^2 * ({p_km1}) + (1)^2 * ({p_kp1}) + (2)^2 * ({p_kp2})")
    print(f"E[Jump^2] for large k = {E_sq_jump_k}")
    print(f"The asymptotic variance is sigma^2 = lim_{{k->oo}} E[Jump^2] = {sigma_sq}\n")
    
    print("Step 3: Apply the Lamperti-Kesten criterion for transience.")
    print("The criterion states the chain is transient if lim_{k->oo} (k * mu_k / sigma^2) > 1/2.\n")

    # Calculate the limit in the criterion
    k_mu_k = sympy.simplify(k * mu_k)
    limit_val = k_mu_k / sigma_sq
    
    print(f"Calculating the terms for the criterion:")
    print(f"k * mu_k = k * ({mu_k}) = {k_mu_k}")
    print(f"limit_val = ({k_mu_k}) / ({sigma_sq}) = {limit_val}\n")

    print(f"Step 4: Solve the inequality {limit_val} > 1/2 for c to find the condition for transience.")

    # Solve the equality to find the boundary/critical value of c
    critical_c_eq = sympy.Eq(limit_val, sympy.Rational(1, 2))
    critical_c_sol = sympy.solve(critical_c_eq, c)
    
    if not critical_c_sol:
        print("Could not solve for a critical value of c.")
        return
        
    critical_c = critical_c_sol[0]

    # Express the final equation with its components
    c_val = sympy.sympify(critical_c)
    num, den = c_val.p, c_val.q

    print(f"Solving the equality {limit_val} = 1/2 gives the boundary:")
    print(f"{4*c}/5 = 1/2")
    print(f"c = {num}/{den}")

    print("\nThe chain is transient for c > 5/8.")
    print("The set of c for which the chain is transient is (5/8, infinity).")
    print(f"Therefore, the infimum of this set is {critical_c}.")


if __name__ == '__main__':
    solve_markov_chain_transience()