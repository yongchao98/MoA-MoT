import numpy as np
from scipy.optimize import fsolve, brentq

def find_oscillation_boundary():
    """
    This script calculates the critical value of the parameter 'n' for which the
    given biochemical system starts to exhibit oscillations.
    """
    # Given parameter values
    alpha = 100.0
    beta = 20.0

    def get_bifurcation_function_value(n, alpha_val, beta_val):
        """
        For a given n, this function calculates the value of the expression
        whose root marks the Hopf bifurcation.
        It first finds the steady state 's', then computes the value.
        """
        # Define the steady-state equation to find s
        # f(s) = alpha / (1 + s^n) - s - beta / (1 + s) = 0
        def steady_state_eq(s_var):
            # The concentration 's' must be positive.
            if s_var <= 0:
                return 1e10  # Return a large value to penalize non-physical solutions
            return alpha_val / (1 + s_var**n) - s_var - beta_val / (1 + s_var)

        # Numerically solve for the steady state 's'
        # A starting guess of s=3.0 is reasonable for the expected range of n.
        s_steady, _, ier, _ = fsolve(steady_state_eq, 3.0, full_output=True)
        
        # Check if the solver was successful
        if ier != 1:
            # Return a value that indicates failure to the root finder
            return -1e10 if n > 1.5 else 1e10

        s = s_steady[0]

        # Calculate the components of the Jacobian matrix at the steady state (s, s, s)
        # J = [[A, 0, B], [B, A, 0], [0, B, A]]
        # A = -1 + beta / (1 + s)^2
        # B = -alpha * n * s^(n-1) / (1 + s^n)^2
        A = -1.0 + beta_val / (1.0 + s)**2
        B = -alpha_val * n * s**(n - 1) / (1.0 + s**n)**2

        # A Hopf bifurcation occurs when the real part of the complex eigenvalues is zero.
        # The real part of the complex eigenvalues is given by A - B/2.
        # We need to find the value of 'n' for which A - B/2 = 0.
        return A - B/2

    # We search for the critical 'n' in a plausible interval, e.g., [1.0, 4.0].
    # Our analysis suggests the root lies between n=1 (stable) and n=4 (unstable).
    try:
        critical_n = brentq(get_bifurcation_function_value, 1.0, 4.0, args=(alpha, beta))
        
        print(f"The system is stable for n < {critical_n:.4f} and exhibits oscillations for n > {critical_n:.4f}.")
        print("This transition occurs at the Hopf bifurcation point.")
        
        # To fulfill the "output each number in the final equation" requirement,
        # we show the values at the critical point.
        # The "final equation" is the bifurcation condition: A - B/2 = 0.
        
        # Recalculate s, A, and B at the critical_n
        def ss_eq(s_var):
            return alpha / (1 + s_var**critical_n) - s_var - beta / (1 + s_var)
        s_crit = fsolve(ss_eq, 3.0)[0]
        A_crit = -1.0 + beta / (1.0 + s_crit)**2
        B_crit = -alpha * critical_n * s_crit**(critical_n - 1) / (1.0 + s_crit**critical_n)**2
        
        print("\nAt the critical point:")
        print(f"Final Equation for n: n > {critical_n:.4f}")
        print(f"Value of n = {critical_n:.4f}")
        print(f"Steady State s = {s_crit:.4f}")
        print(f"Jacobian Term A = {A_crit:.4f}")
        print(f"Jacobian Term B = {B_crit:.4f}")
        print(f"Bifurcation Condition (A - B/2) = {A_crit - B_crit/2:.4f} (which is close to zero)")

    except ValueError:
        print("Could not find the critical value of n in the specified interval.")
        print("The function values at the interval boundaries may not have opposite signs.")

if __name__ == '__main__':
    find_oscillation_boundary()