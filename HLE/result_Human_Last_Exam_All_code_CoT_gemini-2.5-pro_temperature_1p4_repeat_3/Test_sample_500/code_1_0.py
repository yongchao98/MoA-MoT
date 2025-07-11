import sympy

def solve_queueing_problem():
    """
    Solves the queueing problem by analyzing the mean service time.
    """
    # Define the symbolic variable for integration
    u = sympy.Symbol('u', positive=True)
    
    # The arrival rate given in the problem
    lambda_rate = 3
    
    # We are given the tail probability P(S > u) for large u.
    # P(S > u) = 1/(3u) + m/(u*ln(u))
    # We can analyze the convergence of the integral for E[S] by examining its terms.
    # Let's use the first term. If its integral diverges, the whole E[S] integral diverges.
    term_to_analyze = 1 / (3 * u)
    
    # E[S] = integral from 0 to infinity of P(S > u) du.
    # This integral can be split at a large value 'A':
    # E[S] = integral from 0 to A of P(S > u) du + integral from A to infinity of P(S > u) du.
    # The first part is a finite constant. We analyze the second part for convergence.
    A = sympy.Symbol('A', positive=True)
    
    # Use sympy to calculate the integral of the term from A to infinity.
    integral_result = sympy.integrate(term_to_analyze, (u, A, sympy.oo))
    
    print(f"The arrival rate is lambda = {lambda_rate}.")
    print(f"The mean service time E[S] is given by the integral of the tail probability P(S > u).")
    print(f"For large u, P(S > u) contains the term: {term_to_analyze}.")
    print(f"We test the convergence by integrating this term from a large constant 'A' to infinity.")
    
    # The result of the integration will show divergence.
    print(f"\nResult of the integral of {term_to_analyze} from A to infinity is: {integral_result}.")
    
    print("\nSince the integral diverges, the mean service time E[S] is infinite.")
    print("The expected number of customers L in the system is given by lambda * E[S].")
    
    # Final equation showing the numbers involved.
    print(f"L = {lambda_rate} * E[S] = {lambda_rate} * infinity = infinity.")
    
    print("\nBecause the expected number of customers is infinite, the number of customers in the system, X_t, grows without bound as t -> infinity.")
    print("Therefore, the limit inferior of X_t is also infinity.")

solve_queueing_problem()
