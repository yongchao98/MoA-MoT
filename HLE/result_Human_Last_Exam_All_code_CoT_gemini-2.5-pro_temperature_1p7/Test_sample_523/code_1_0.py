import numpy as np
from scipy.optimize import root

def find_oscillation_boundary():
    """
    Calculates the critical value of the parameter 'n' at which the system
    starts to exhibit oscillations.
    """
    # Parameters of the system from the problem description
    alpha = 100.0
    beta = 20.0

    def hopf_bifurcation_eqs(variables):
        """
        Defines the system of two non-linear equations that must be solved
        to find the Hopf bifurcation point.
        
        Args:
            variables (list): A list containing [s, n], where 's' is the
                              steady-state concentration and 'n' is the Hill coefficient.
        
        Returns:
            list: The values of the two equations, which should be zero at the solution.
        """
        s, n = variables
        
        # Ensure s and n are positive to avoid math errors in the equations
        if s <= 0 or n <= 0:
            return [1e9, 1e9] # Return large error for invalid domain

        try:
            # Equation 1: The steady-state condition
            # 0 = alpha / (1 + s^n) - s - beta / (1 + s)
            common_term_s_n = 1.0 + s**n
            eq1 = alpha / common_term_s_n - s - beta / (1.0 + s)

            # Equation 2: The Hopf bifurcation condition (Real part of eigenvalue is zero)
            # Re(lambda) = a - b/2 = 0
            # where a = -1 + beta / (1+s)^2
            # and b = -alpha * n * s^(n-1) / (1+s^n)^2
            a = -1.0 + beta / ((1.0 + s)**2)
            b = -alpha * n * s**(n - 1.0) / (common_term_s_n**2)
            eq2 = a - b / 2.0
            
            # Check for non-real results which can happen with invalid inputs
            if isinstance(eq1, complex) or isinstance(eq2, complex):
                 return [1e9, 1e9]
                 
        except (ValueError, OverflowError, ZeroDivisionError):
            return [1e9, 1e9] # Return large error if calculation fails

        return [eq1, eq2]

    # Provide an initial guess for the solver [s_guess, n_guess]
    initial_guess = [4.0, 1.8]

    # Use a numerical root finder to solve the system of equations
    solution = root(hopf_bifurcation_eqs, initial_guess, method='hybr')

    if solution.success:
        s_crit, n_crit = solution.x
        print("The analysis of the system's stability leads to the following conclusion:")
        print(f"The system is stable for low values of n and transitions to oscillations via a Hopf bifurcation.")
        print(f"The critical point for this transition occurs at n ≈ {n_crit:.3f}.")
        print("\nTherefore, the final equation describing the condition for oscillations is:")
        print(f"n > {n_crit:.3f}")
    else:
        print("The numerical solver failed to find the critical point.")
        print(f"Solver message: {solution.message}")

if __name__ == '__main__':
    find_oscillation_boundary()