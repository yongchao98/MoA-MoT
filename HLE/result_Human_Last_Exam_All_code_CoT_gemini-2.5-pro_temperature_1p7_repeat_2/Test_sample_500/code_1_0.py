import sympy

def analyze_queueing_system():
    """
    Analyzes the described M/G/infinity queueing system to find the
    limit inferior of the number of customers in the system as t -> infinity.
    """

    print("Step 1: Identifying the Queueing System Model")
    print("---------------------------------------------")
    print("Arrivals: Poisson process (rate lambda=3). This corresponds to 'M' (Markovian).")
    print("Service Times: General i.i.d. distribution. This corresponds to 'G' (General).")
    print("Servers: Customers enter service immediately, implying no queue. This means there is always a free server, corresponding to 'infinity' servers.")
    print("Conclusion: The system is an M/G/∞ queue.\n")

    print("Step 2: Condition for Stability")
    print("-------------------------------")
    print("The stability of an M/G/∞ queue depends on the expected service time, E[S].")
    print("If E[S] is finite, the system is stable.")
    print("If E[S] is infinite, the number of customers in the system grows without bound.\n")

    print("Step 3: Calculating Expected Service Time E[S]")
    print("---------------------------------------------")
    print("The expected service time E[S] is the integral of the survival function P(S > u) from 0 to infinity.")
    print("E[S] = integral from 0 to infinity of P(S > u) du.")
    print("We are given that for large enough u, P(S > u) = 1/(3*u) + m/(u*ln(u)).")
    print("To determine if E[S] is finite, we only need to check if the integral of this tail converges.\n")

    print("Step 4: Using SymPy to Analyze the Integral's Convergence")
    print("---------------------------------------------------------")
    # Define symbolic variables
    u = sympy.Symbol('u', positive=True)
    A = sympy.Symbol('A', positive=True) # A represents a "large enough" number

    # The survival function has two terms. We only need to show one of them diverges.
    # The term 1/(3*u) is simpler and sufficient.
    integrand_term = 1 / (3 * u)

    print(f"We will check the convergence of the integral of the term: {integrand_term}")
    
    # Calculate the definite integral from A to infinity
    integral_result = sympy.integrate(integrand_term, (u, A, sympy.oo))
    
    print(f"The integral of {integrand_term} from A to infinity is: {integral_result}")

    if integral_result == sympy.oo:
        print("The integral diverges. This means that E[S] = ∞.\n")
    else:
        print("The integral converges. This would mean E[S] is finite.\n")

    print("Step 5: Deducing the System's Long-Term Behavior")
    print("--------------------------------------------------")
    print("Since E[S] is infinite, the M/G/∞ system is unstable.")
    print("This implies that the number of customers in the system, X_t, grows without bound as t -> infinity.")
    print("In mathematical terms, lim (t->∞) X_t = ∞ (almost surely).\n")

    print("Step 6: Final Calculation of liminf")
    print("------------------------------------")
    print("The limit inferior (liminf) of a quantity that tends to infinity is also infinity.")
    final_answer = sympy.oo
    print(f"Therefore, liminf (t->∞) X_t = {final_answer}")


if __name__ == '__main__':
    analyze_queueing_system()