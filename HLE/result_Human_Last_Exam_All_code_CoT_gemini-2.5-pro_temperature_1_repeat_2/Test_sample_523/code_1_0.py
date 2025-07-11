import numpy as np
from scipy.optimize import fsolve

def find_oscillatory_range():
    """
    Determines the critical value of n for the onset of oscillations
    in the given biochemical system.
    """
    # Given parameters
    alpha = 100
    beta = 20

    print(f"Analyzing the system with parameters alpha = {alpha} and beta = {beta}")

    def find_steady_state(n, alpha_val, beta_val):
        """Numerically finds the steady state concentration 's' for a given n."""
        # We need to solve the equation: alpha / (1 + s^n) - s - beta / (1 + s) = 0
        steady_state_eq = lambda s: alpha_val / (1 + s**n) - s - beta_val / (1 + s)
        
        # An initial guess of s=5.0 is robust for the expected range of n.
        s_initial_guess = 5.0
        s_steady, info, ier, msg = fsolve(steady_state_eq, s_initial_guess, full_output=True)
        
        # Ensure the root finding was successful
        if ier != 1:
            raise RuntimeError(f"Failed to find a steady state for n={n}. Message: {msg}")
        
        return s_steady[0]

    def hopf_bifurcation_condition(n, alpha_val, beta_val):
        """
        Calculates the real part of the complex conjugate eigenvalues.
        A value of 0 indicates the Hopf bifurcation point.
        """
        # The parameter n must be positive.
        if n <= 0:
            return -1  # Return a value indicating stability

        try:
            # First, find the steady state 's' for the given 'n'
            s = find_steady_state(n, alpha_val, beta_val)
            if s <= 0:
                return -1 # Unphysical steady state
        except RuntimeError:
            return -1 # Failed to find a solution, assume stability

        # Calculate the terms of the Jacobian matrix evaluated at the steady state
        A = -1 + beta_val / (1 + s)**2
        B = -alpha_val * n * s**(n - 1) / (1 + s**n)**2
        
        # The real part of the complex eigenvalues is A - B/2.
        # Oscillations occur when this term is > 0.
        return A - B/2

    # We need to find the root of the Hopf condition function to find the critical n.
    # From preliminary analysis, the root is expected to be slightly above 1.
    initial_guess_for_n = 1.2
    critical_n, info, ier, msg = fsolve(
        hopf_bifurcation_condition, 
        initial_guess_for_n, 
        args=(alpha, beta),
        full_output=True
    )

    if ier != 1:
        print("\nCould not reliably determine the critical value of n.")
        return

    n_c = critical_n[0]

    print("\nThe stability of the system's steady state depends on the parameter n.")
    print("Oscillations occur when the steady state becomes unstable via a Hopf bifurcation.")
    print(f"We have numerically solved for the critical value of n where this bifurcation occurs.")
    print("\nFinal Equation for Oscillation Threshold:")
    print(f"For oscillations to occur with alpha = {alpha} and beta = {beta}, the condition is:")
    print(f"n > {n_c:.4f}")
    
    # Return the final answer in the required format
    print(f"\n<<<{n_c:.4f}>>>")

find_oscillatory_range()