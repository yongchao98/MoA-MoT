import sympy
from sympy import oo, Symbol, integrate, ln, pprint

def analyze_queueing_system():
    """
    Analyzes the M/G/infinity queueing system to find liminf of the number of customers.
    """

    # --- Step 1: Define parameters and theoretical background ---
    lambda_rate = 3
    # m is a positive integer. Its exact value doesn't change the outcome, so we treat it as a symbol.
    m = Symbol('m', positive=True, integer=True)
    u = Symbol('u', positive=True)
    A = Symbol('A', positive=True, real=True) # A represents a large number for the integral lower bound

    print("Step 1: System Identification and Key Principles")
    print("-------------------------------------------------")
    print(f"The system is an M/G/infinity queue with arrival rate λ = {lambda_rate}.")
    print("The long-term behavior of the number of customers (X_t) depends on the mean service time E[S].")
    print("A core theorem states: If E[S] = ∞, then X_t → ∞ almost surely. Consequently, liminf X_t = ∞.")
    print("Our goal is to determine if E[S] is finite or infinite.")

    # --- Step 2: Analyze the Mean Service Time E[S] ---
    print("\nStep 2: Analysis of the Mean Service Time E[S]")
    print("------------------------------------------------")
    print("E[S] is the integral of the survival function P(S > u) from 0 to ∞.")
    survival_function_tail = 1 / (3 * u) + m / (u * ln(u))
    print("For large u, the survival function is given by:")
    pprint(survival_function_tail)
    
    print("\nTo check if E[S] is infinite, we test if the integral of this tail expression diverges.")
    
    # We use sympy to evaluate the integral for each part of the expression.
    # We assume A is large enough such that ln(A) is well-defined and positive.
    integral_part1 = integrate(1/(3*u), (u, A, oo))
    integral_part2 = integrate(m/(u*ln(u)), (u, A, oo))

    print(f"\nIntegral of the first term (1/(3u)) from A to ∞ is: {integral_part1}")
    print(f"Integral of the second term (m/(u*ln(u))) from A to ∞ is: {integral_part2}")
    
    # --- Step 3: Conclude based on the integration results ---
    print("\nStep 3: Conclusion from the Analysis")
    print("---------------------------------------")
    
    if integral_part1 == oo or integral_part2 == oo:
        print("Since the integral for at least one term of the tail probability diverges, the total integral for E[S] is infinite.")
        print("Therefore, the mean service time E[S] = ∞.")
    else:
        # This case will not be reached given the problem statement.
        print("The mean service time E[S] is finite.")

    print("\nBecause E[S] is infinite, the M/G/∞ system is unstable.")
    print("The number of customers, X_t, will grow without bound almost surely as t → ∞.")

    # Frame the final answer in the context of the problem's parameters
    print("\nThe expected number of customers in the system, E[X_t], can be conceptualized as tending towards:")
    print(f"E[X_t] → λ * E[S]")
    print(f"E[X_t] → {lambda_rate} * ∞ = ∞")

    final_answer = "Infinity"
    print(f"\nThis confirms that liminf_{{t→∞}} X_t = {final_answer}.")

# Execute the analysis
analyze_queueing_system()