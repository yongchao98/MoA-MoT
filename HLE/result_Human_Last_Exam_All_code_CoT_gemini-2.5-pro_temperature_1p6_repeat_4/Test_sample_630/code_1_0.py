import numpy as np

def calculate_errors(T, D):
    """
    Calculates optimization, approximation, and total error based on a trade-off.

    This function simulates the trade-off analysis for stochastic logistic regression.
    It finds the optimal internal radius R to trade-off optimization error
    and approximation error, and then computes the resulting total error.
    """

    # We model the error components as functions of an internal radius R.
    # OptErr(R) = c1 * exp(R) / T
    # ApproxErr(R) = c2 * exp(-R)
    # We set constants c1 and c2 to 1 for simplicity.

    # Find the optimal radius R by minimizing OptErr + ApproxErr.
    # The minimum occurs when the two terms are equal: exp(R)/T = exp(-R) => R = 0.5*log(T)
    R_optimal = 0.5 * np.log(T)

    if R_optimal > D:
        # If the theoretically optimal R is outside our main domain, we are limited by D.
        # This can happen if T is very large compared to exp(D).
        # In the regime T = O(e^D), R_optimal <= D holds for large T.
        print(f"Warning: Optimal radius R={R_optimal:.2f} is larger than D={D:.2f}. Capping at R=D.")
        R_optimal = D

    # Calculate the error components at the optimal R.
    opt_error = np.exp(R_optimal) / T
    # The approximation error model L(w_R) - L(w_D) simplifies to exp(-R) for large D.
    approx_error = np.exp(-R_optimal)

    total_error = opt_error + approx_error
    
    # The theoretical optimal rate is Theta(1/sqrt(T))
    theoretical_rate = 1.0 / np.sqrt(T)

    print(f"For T = {T}:")
    print(f"Domain radius D is set to {D:.2f} (since T is in the order of e^D).")
    print(f"Optimal internal radius R to trade-off errors is {R_optimal:.2f}.")
    print(f"Approximation Error at R_optimal: {approx_error:.6f}")
    print(f"Optimization Error at R_optimal:   {opt_error:.6f}")
    print(f"Total Minimized Error:              {total_error:.6f}")
    print(f"Theoretical rate O(1/sqrt(T)):      {theoretical_rate:.6f}")
    print("-" * 20)


# Run simulation for different values of T
T_values = [1e4, 1e6, 1e8]
for T in T_values:
    # Set D such that T is in the order of e^D, i.e., D is in the order of log(T)
    D = np.log(T)
    calculate_errors(T, D)
