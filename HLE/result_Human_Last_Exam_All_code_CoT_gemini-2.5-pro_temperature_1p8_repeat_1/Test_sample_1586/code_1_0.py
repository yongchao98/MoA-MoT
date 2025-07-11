import sympy as sp

def solve_markov_problem():
    """
    Solves for the supremum of alpha for which the alpha-th moment of the hitting time is finite.
    """
    # Define symbolic variables for n, c, and alpha.
    # n is the state, c is the parameter from the problem, alpha is the moment exponent.
    n, c, alpha = sp.symbols('n c alpha', positive=True, real=True)

    print("Step 1: Define the transition probability ratio rho_n for large n.")
    # For large n, p(n, n+1) = 1/2 - c/n
    # And p(n, n-1) = 1 - p(n, n+1) = 1/2 + c/n
    p_plus = sp.Rational(1, 2) - c / n
    p_minus = sp.Rational(1, 2) + c / n
    
    # rho_n is the ratio p(n, n-1) / p(n, n+1)
    rho_n = p_minus / p_plus
    rho_n_simplified = sp.simplify(rho_n)
    print(f"rho_n = {rho_n_simplified}\n")

    print("Step 2: Find the asymptotic behavior of rho_n for large n.")
    # The series expansion of rho_n around n = infinity shows its behavior.
    rho_n_series = sp.series(rho_n_simplified, n, sp.oo, 2)
    print(f"For large n, rho_n is approximately: {rho_n_series}\n")

    print("Step 3: Find the asymptotic behavior of the product term in the criterion.")
    # We analyze the product by looking at the sum of the logarithms.
    # log(rho_n) for large n is approximately:
    log_rho_n_series = sp.series(sp.log(rho_n_simplified), n, sp.oo, 2)
    leading_term_log_rho = log_rho_n_series.args[0]
    print(f"For large n, log(rho_n) is approximately: {leading_term_log_rho}")
    
    # The sum Sum[log(rho_k)] from k=1 to n behaves like the integral of the leading term.
    # The integral of (4*c/n) is 4*c*log(n).
    # So, log(Product[rho_k]) is asymptotic to 4*c*log(n).
    # This implies Product[rho_k] is asymptotic to exp(4*c*log(n)) = n^(4*c).
    asymptotic_product = n**(4 * c)
    print(f"The product Product[rho_k for k=1..n] is asymptotic to: {asymptotic_product}\n")
    
    print("Step 4: Apply the convergence criterion.")
    # The criterion requires the sum of n^(alpha-1) * (product)^(-1) to converge.
    # The general term of the series is therefore asymptotic to:
    series_term = n**(alpha - 1) * asymptotic_product**(-1)
    series_term = sp.simplify(series_term)
    print(f"The general term of the criterion series is asymptotic to: {series_term}")
    
    # This is a p-series of the form Sum[n**(-p)]. It converges if p > 1.
    p = -1 * (alpha - 1 - 4 * c)
    print(f"This is a p-series with p = -(alpha - 1 - 4*c) = {p}")
    
    # The condition for convergence is p > 1.
    convergence_condition = sp.Gt(p, 1)
    print(f"The convergence condition is p > 1, which means: {convergence_condition}\n")

    print("Step 5: Solve for alpha to find the supremum.")
    # We solve the inequality for alpha.
    alpha_solution = sp.solve(convergence_condition, alpha)
    print(f"The condition for E[tau^alpha] to be finite is: {alpha_solution}")
    
    # The supremum is the boundary of this set.
    final_answer = alpha_solution.rhs
    
    print("\n--- FINAL RESULT ---")
    print(f"The supremum of alpha is the expression: {final_answer}")
    
    # As requested, output the components of the final expression.
    coeff, param = final_answer.as_coeff_mul()
    print(f"In the final equation sup(alpha) = {final_answer}:")
    print(f"The number in the expression is: {coeff}")
    print(f"The parameter in the expression is: '{param[0]}'")

solve_markov_problem()