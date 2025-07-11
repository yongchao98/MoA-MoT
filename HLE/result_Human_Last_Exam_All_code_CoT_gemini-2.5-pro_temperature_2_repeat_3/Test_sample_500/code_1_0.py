import numpy as np
import scipy.integrate as spi
import warnings

def solve_queueing_problem():
    """
    This script solves the given queueing theory problem by analyzing the
    system's properties and calculating the limit inferior of the number of customers.
    """

    # --- Problem Parameters ---
    lambda_rate = 3
    # The problem states m is a positive integer. Its exact value (as long as > 0)
    # does not affect the divergence, so we use m=1 for the analysis.
    m = 1

    print("--- Step 1: System Identification ---")
    print("The system is an M/G/infinity queue.")
    print(" - Arrivals are a Poisson process ('M') with rate lambda.")
    print(" - Service times have a general distribution ('G').")
    print(" - Customers enter service immediately, implying infinite servers ('infinity').")
    print(f"The given arrival rate is lambda = {lambda_rate}.\n")

    print("--- Step 2: Analyze Expected Service Time E[S] ---")
    print("The stability of an M/G/infinity queue depends on the expected service time, E[S].")
    print("E[S] is calculated by integrating the tail probability P(S > u) from 0 to infinity.")
    print(f"For large u, the problem states: P(S > u) = 1/(3*u) + {m}/(u*ln(u)).")
    print("We must check if the integral of this function from 0 to infinity converges or diverges.")

    # We can determine convergence by analyzing the integral of the terms for large u.
    # Analytically, the integral of m/(u*ln(u)) is m*ln(ln(u)).
    # As u -> infinity, ln(ln(u)) also -> infinity. This means the integral diverges.
    # Therefore, the expected service time E[S] is infinite.
    # Let's numerically demonstrate this divergence to support the analytical conclusion.

    def p_s_gt_u(u, m_val=1):
        # The formula is for "large enough u". Let's assume this holds for u > 1.
        if u <= 1.000001:  # Avoid u=1 where ln(u)=0
            return 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return (1 / (3 * u)) + (m_val / (u * np.log(u)))

    print("\nDemonstrating the divergence of the integral part of E[S]:")
    lower_bound = 2.0
    upper_bounds = [1e3, 1e6, 1e9]
    for T in upper_bounds:
        integral_val, _ = spi.quad(lambda u: p_s_gt_u(u, m), lower_bound, T)
        print(f"  - Integral of P(S > u) from {lower_bound} to {T:.0e}: {integral_val:.4f}")

    print("\nAs the integration's upper bound increases, the resulting value grows without limit.")
    print("This confirms that the integral diverges, and therefore, E[S] = infinity.\n")

    print("--- Step 3: Determine Long-Term System Behavior ---")
    print("For an M/G/infinity queue, the expected number of customers in the system, E[X_t],")
    print("tends towards lambda * E[S] as t -> infinity.")
    print(f"In this case, E[X_t] tends to {lambda_rate} * infinity, which is infinity.")
    print("Since the expected number of customers grows without bound, the actual number of customers, X_t, must also tend to infinity.\n")

    print("--- Step 4: Final Conclusion ---")
    print("The limit inferior, liminf as t->infinity of X_t, is the smallest value that X_t gets arbitrarily close to for all large t.")
    print("Since X_t grows to infinity, it will eventually exceed any finite number and never return below it.")
    print("Therefore, its limit inferior must be infinity.\n")
    
    print("Final logical flow:")
    print(f"  Arrival Rate (lambda) = {lambda_rate}")
    print(f"  Service Time Tail Probability (P(S > u)) for large u contains a term {m}/(u*ln(u))")
    print(f"  Expected Service Time (E[S]) = Integral from 0 to inf of P(S > u) du, which diverges. So, E[S] = infinity")
    print(f"  Expected number of customers (E[X_t]) = lambda * E[S] = {lambda_rate} * infinity. So, E[X_t] -> infinity")
    print(f"  Number of customers (X_t) -> infinity")


if __name__ == '__main__':
    solve_queueing_problem()