import sympy
from sympy import oo, log, integrate, Symbol

def solve_queueing_problem():
    """
    This script solves the queueing problem by analyzing the expected service time.
    """
    # Step 1: Define the system parameters from the problem description.
    lambda_rate = 3
    u = Symbol('u', positive=True)
    m = Symbol('m', positive=True, integer=True)

    print("Step 1: Identifying the Queueing Model")
    print("The system described is an M/G/infinity queue.")
    print("- Arrivals follow a Poisson process (M) with rate lambda = 3.")
    print("- Service times have a general distribution (G).")
    print("- Customers enter service immediately, implying an infinite number of servers (infinity).")
    print("-" * 50)

    print("Step 2: Determining the Condition for Long-Term Behavior")
    print("In an M/G/infinity queue, the number of customers in the system, X_t, grows to infinity if the expected service time, E[S], is infinite.")
    print("We need to calculate E[S] to determine the behavior of X_t as t -> infinity.")
    print("-" * 50)

    print("Step 3: Calculating the Expected Service Time E[S]")
    print("The expected service time E[S] is the integral of the tail probability P(S > u) from 0 to infinity.")
    print("E[S] = integral(P(S > u) du) from 0 to oo")
    print("For large u, the tail probability is given as:")
    tail_prob_formula = 1/(3*u) + m/(u * log(u))
    print(f"P(S > u) = {tail_prob_formula}")
    print("To check if E[S] is finite, we only need to analyze the convergence of the integral of this expression to infinity.")
    print("-" * 50)

    print("Step 4: Evaluating the Integral for Divergence")
    # We analyze the two parts of the integral separately.
    # The integral of 1/(3*u) is (1/3)*log(u).
    # The integral of m/(u*log(u)) is m*log(log(u)).
    # Both terms diverge as u -> infinity.
    
    # Let's use sympy to show this formally.
    antiderivative = integrate(tail_prob_formula, u)
    print(f"The indefinite integral of the tail probability is: {antiderivative}")
    
    limit_at_infinity = sympy.limit(antiderivative, u, oo)
    print(f"The limit of this antiderivative as u -> infinity is: {limit_at_infinity}")
    print("Since the integral diverges, the expected service time E[S] is infinite.")
    print("-" * 50)

    print("Step 5: Final Conclusion")
    print("Because E[S] is infinite, the number of customers in the system, X_t, grows without bound almost surely.")
    print("This means that for any large number K, there is a time T such that for all t > T, X_t > K.")
    print("Therefore, the limit inferior of X_t as t approaches infinity is also infinity.")
    
    # The final equation is liminf_{t->inf} X_t = infinity.
    # The number from the problem is lambda = 3.
    print("\nFinal Equation:")
    print(f"Given lambda = {lambda_rate}, the final result is:")
    print("liminf_{t->oo} X_t = infinity")

solve_queueing_problem()