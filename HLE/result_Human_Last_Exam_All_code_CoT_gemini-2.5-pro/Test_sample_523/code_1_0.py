import numpy as np
from scipy.optimize import fsolve, brentq

def find_critical_n():
    """
    This script finds the critical value of the parameter 'n' for which the
    given biochemical system starts to exhibit oscillations.
    """
    # Define system parameters from the problem description
    alpha = 100.0
    beta = 20.0

    def get_stability_metric(n, alpha_val, beta_val):
        """
        Calculates the stability metric (real part of the critical eigenvalue) for a given n.
        The system becomes unstable (oscillatory) when this metric crosses from negative to positive.
        A positive value indicates oscillations.
        """
        # Define the steady-state equation to be solved for x_s
        def steady_state_eq(x):
            # We need x > 0 for the terms to be real.
            if x <= 0:
                return float('inf')
            return alpha_val / (1 + x**n) - x - beta_val / (1 + x)

        # Numerically solve for the steady-state concentration x_s
        # We use an initial guess of x=2.0, as it's in a reasonable range.
        try:
            x_s, = fsolve(steady_state_eq, 2.0, xtol=1e-9)
        except Exception:
            # If the solver fails, it's not a valid physical state.
            return -1 # Return a value indicating stability

        # If the solver returns a non-physical value, treat as stable.
        if x_s <= 0:
            return -1

        # Calculate the Jacobian elements 'a' and 'b' at the steady state
        a = -1.0 + beta_val / (1 + x_s)**2
        b = - (alpha_val * n * x_s**(n - 1)) / (1 + x_s**n)**2
        
        # The stability is determined by the sign of 'a - b/2'.
        # This is the real part of the pair of complex eigenvalues that lead to oscillation.
        return a - b / 2.0

    # We need to find the root of `get_stability_metric(n) = 0`.
    # Based on preliminary analysis, the root lies between n=1 and n=2.
    # `brentq` is a robust root-finding algorithm for a bracketed root.
    try:
        critical_n = brentq(get_stability_metric, 1.0, 2.0, args=(alpha, beta))
    except ValueError:
        print("Error: Could not find the critical value of n in the specified interval [1, 2].")
        print("The stability metric may not cross zero in this range.")
        return

    # Now, calculate the values at the critical point for the final output
    x_s_crit, = fsolve(lambda x: alpha / (1 + x**critical_n) - x - beta / (1 + x), 2.0)
    a_crit = -1.0 + beta / (1 + x_s_crit)**2
    b_crit = - (alpha * critical_n * x_s_crit**(critical_n - 1)) / (1 + x_s_crit**critical_n)**2

    # Print the results
    print("The system exhibits oscillations for values of n greater than a critical value, n_c.")
    print(f"The critical value n_c where oscillations begin is: {critical_n:.4f}")
    print("\nThis critical value is found by solving the bifurcation condition equation.")
    print("The final equation for the bifurcation point is: a - b/2 = 0")
    print("\nAt the point of bifurcation:")
    print(f"The steady-state concentration is x_s = {x_s_crit:.4f}")
    print(f"The value of the Jacobian element 'a' is {a_crit:.4f}")
    print(f"The value of the Jacobian element 'b' is {b_crit:.4f}")
    print("\nSubstituting these numbers into the final equation:")
    print(f"a - b/2 = {a_crit:.4f} - ({b_crit:.4f}) / 2 = {a_crit - b_crit/2:.6f}")
    print("This value is approximately zero, confirming our result.")

if __name__ == '__main__':
    find_critical_n()
<<<1.2859>>>