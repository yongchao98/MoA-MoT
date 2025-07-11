import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    """
    # Define symbols for state k and parameter c
    k = sympy.Symbol('k')
    c = sympy.Symbol('c')

    # Define transition probabilities for large k as per the problem description
    # P(k -> k+1)
    p_plus_1 = sympy.Rational(1, 4) + c/k
    # P(k -> k-1)
    p_minus_1 = sympy.Rational(1, 4) - c/k
    # P(k -> k+2)
    p_plus_2 = sympy.Rational(1, 4)
    # P(k -> k-2)
    p_minus_2 = sympy.Rational(1, 4)

    # Step 1: Calculate the drift (mean displacement) mu_k
    # mu_k = E[X_{n+1} - X_n | X_n=k]
    # The displacements are +1, -1, +2, -2
    mu_k = (1 * p_plus_1) + (-1 * p_minus_1) + (2 * p_plus_2) + (-2 * p_minus_2)
    mu_k_simplified = sympy.simplify(mu_k)

    print("Step 1: Calculate the drift (mean displacement) mu_k.")
    print(f"mu_k = (1) * (1/4 + c/k) + (-1) * (1/4 - c/k) + (2) * (1/4) + (-2) * (1/4)")
    print(f"Simplified mu_k = {mu_k_simplified}\n")

    # Step 2: Calculate the second moment of displacement sigma_sq_k
    # sigma_sq_k = E[(X_{n+1} - X_n)^2 | X_n=k]
    sigma_sq_k = ((1)**2 * p_plus_1) + ((-1)**2 * p_minus_1) + ((2)**2 * p_plus_2) + ((-2)**2 * p_minus_2)
    sigma_sq_k_simplified = sympy.simplify(sigma_sq_k)

    print("Step 2: Calculate the second moment of displacement sigma_sq_k.")
    print(f"sigma_sq_k = (1)^2 * (1/4 + c/k) + (-1)^2 * (1/4 - c/k) + (2)^2 * (1/4) + (-2)^2 * (1/4)")
    print(f"Simplified sigma_sq_k = {sigma_sq_k_simplified}\n")

    # Step 3: Apply the Pakes-Lamperti criterion for transience.
    # The chain is transient if alpha = lim_{k->inf} (k * 2 * mu_k / sigma_sq_k) > 1.
    alpha_expr = k * 2 * mu_k_simplified / sigma_sq_k_simplified
    alpha = sympy.limit(alpha_expr, k, sympy.oo)
    
    print("Step 3: Apply the transience criterion.")
    print("The chain is transient if alpha = lim_{k->inf} (k * 2 * mu_k / sigma_sq_k) > 1.")
    print(f"The expression for alpha is: {alpha}\n")

    # Step 4: Find the infimum of c by solving alpha = 1.
    # The set of c for which the chain is transient is {c | alpha > 1}.
    # The infimum of this set is the value of c where alpha = 1.
    equation_to_solve = sympy.Eq(alpha, 1)
    solution = sympy.solve(equation_to_solve, c)
    infimum_c = solution[0]
    
    # Print the final equation with numbers as requested
    lhs = str(equation_to_solve.lhs).replace("c", "*c")
    rhs = str(equation_to_solve.rhs)
    
    print("Step 4: Find the infimum of c by solving for the boundary case alpha = 1.")
    print(f"The equation to solve is: {lhs} = {rhs}")
    print(f"The solution for c, which is the required infimum, is: {infimum_c}")

if __name__ == '__main__':
    solve_markov_chain_transience()
    print("\n<<<5/8>>>")