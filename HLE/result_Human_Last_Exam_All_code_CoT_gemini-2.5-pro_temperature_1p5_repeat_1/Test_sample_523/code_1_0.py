import numpy as np
from scipy.optimize import fsolve

def find_bifurcation_point():
    """
    This function solves for the critical values of the fixed point (s) and
    cooperativity (n) at which a Hopf bifurcation occurs in the repressilator system.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    def system_equations(variables):
        """
        Defines the system of two non-linear equations to find the bifurcation point.
        Args:
            variables (list): A list containing the two variables [s, n].
        Returns:
            list: The residuals of the two equations.
        """
        s, n = variables
        
        # Ensure variables are in a valid domain to prevent math errors.
        if s <= 0 or n <= 0:
            return [1e6, 1e6]

        # Equation 1: The fixed-point condition.
        # α / (1 + s^n) - s - β / (1 + s) = 0
        try:
            term_s_n = s**n
        except OverflowError:
            return [1e6, 1e6]
            
        fixed_point_eq = alpha / (1 + term_s_n) - s - beta / (1 + s)

        # Equation 2: The Hopf bifurcation condition.
        # A - B/2 = 0, which is equivalent to 2*A - B = 0
        A = -1.0 + beta / (1 + s)**2
        B = -alpha * n * s**(n - 1) / (1 + term_s_n)**2
        bifurcation_eq = A - B / 2.0

        return [fixed_point_eq, bifurcation_eq]

    # Provide an initial guess for the solver [s_guess, n_guess].
    # Based on preliminary analysis, s is expected to be around 2-4 and n around 2-3.
    initial_guess = [3.0, 2.5]
    
    # Use fsolve to find the solution.
    s_critical, n_critical = fsolve(system_equations, initial_guess)

    # Output the results in a readable format.
    print(f"Given parameters: alpha = {alpha}, beta = {beta}")
    print("\nThe final equation describing the condition for oscillations is derived from the Hopf bifurcation point.")
    print(f"The critical value for the cooperativity 'n' is found to be: {n_critical:.4f}")
    print(f"The system exhibits oscillations for values of n greater than {n_critical:.4f}.")

find_bifurcation_point()