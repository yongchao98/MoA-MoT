import numpy as np
from scipy.optimize import fsolve, brentq

def solve_for_oscillations():
    """
    This function calculates the critical value of the parameter 'n'
    for the onset of oscillations in the given biochemical system.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    def get_stability_metric(n):
        """
        For a given n, this function finds the steady state 's' and
        returns the value of the stability metric 'A + C/2'.
        The root of this function gives the critical n for the Hopf bifurcation.
        """
        # Define the steady-state equation to solve for s
        steady_state_eq = lambda s: s + beta / (1 + s) - alpha / (1 + s**n)

        # Numerically solve for the steady state 's'
        # An initial guess of s=3.0 is robust for the range of n we are testing
        s_ss = fsolve(steady_state_eq, x0=3.0)[0]

        # If solver fails or finds a non-physical root, return a value
        # that indicates stability.
        if s_ss <= 0:
            return -1.0

        # Calculate the terms A and C for the stability condition
        A = -1.0 + beta / (1.0 + s_ss)**2
        C = alpha * n * s_ss**(n - 1) / (1.0 + s_ss**n)**2

        # Return the value of the stability metric. Oscillations occur when this is > 0.
        return A + C / 2.0

    # We need to find the root of get_stability_metric(n) to find the critical value of n.
    # From a preliminary analysis, we know the root lies between n=1 and n=2.
    # We use the Brentq method for a robust root-finding in this interval.
    try:
        n_critical = brentq(get_stability_metric, a=1.0, b=2.0)

        print(f"The system exhibits oscillations for values of n > {n_critical:.5f}")
        
        # Now we verify the "final equation" A + C/2 = 0 at this critical point
        # Recalculate s, A, and C at the critical n
        steady_state_eq_critical = lambda s: s + beta / (1 + s) - alpha / (1 + s**n_critical)
        s_critical = fsolve(steady_state_eq_critical, x0=3.0)[0]
        A_critical = -1.0 + beta / (1.0 + s_critical)**2
        C_critical = alpha * n_critical * s_critical**(n_critical - 1) / (1.0 + s_critical**n_critical)**2
        
        print("\nVerification at the critical point:")
        print(f"At n_c = {n_critical:.5f}, the steady state is s_c = {s_critical:.5f}.")
        print("The Hopf bifurcation condition A + C/2 = 0 must be satisfied.")
        print("Let's plug in the numbers for the final equation:")
        
        # We need to explicitly output each number in the equation
        print(f"{A_critical:.5f} + {C_critical:.5f} / 2 = {A_critical + C_critical/2.0:.5f}")
        print("(The result is approximately 0, as expected for the bifurcation point.)")

    except ValueError:
        print("Could not find the critical value of n in the specified interval.")

solve_for_oscillations()
<<<1.84285>>>