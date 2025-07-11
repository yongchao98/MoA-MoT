import numpy as np
from scipy.optimize import fsolve

def find_critical_n():
    """
    Calculates the critical value of the parameter 'n' for the onset of oscillations
    in the given biochemical system.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    def steady_state_eq(x, n):
        """Equation for the steady state x_ss."""
        return alpha / (1 + x**n) - x - beta / (1 + x)

    def hopf_bifurcation_condition(n):
        """
        This function represents the Hopf bifurcation condition, Re(lambda) = A - B/2 = 0.
        We will find the root of this function, which is the critical value of n.
        """
        # First, find the steady state x_ss for the given n.
        # An initial guess of x0=3.5 is based on preliminary analysis.
        try:
            x_ss = fsolve(steady_state_eq, x0=3.5, args=(n))[0]
            if x_ss <= 0:
                # Return a large value if the steady state is not physically meaningful
                return np.inf
        except (RuntimeError, ValueError):
            # If the solver fails, return a large value
            return np.inf

        # Calculate the terms A and B from the Jacobian matrix
        A = -1 + beta / (1 + x_ss)**2
        B = -alpha * n * x_ss**(n - 1) / (1 + x_ss**n)**2
        
        # The condition for the Hopf bifurcation is A - B/2 = 0
        return A - B / 2

    # Find the value of n that satisfies the Hopf condition using a root-finding algorithm.
    # An initial guess of n=4.0 is used.
    try:
        n_critical = fsolve(hopf_bifurcation_condition, x0=4.0)[0]
    except (RuntimeError, ValueError):
        n_critical = "Error: Could not determine the critical value."

    # Output the parameters and the result
    print(f"For the given system with parameters alpha = {alpha} and beta = {beta},")
    print(f"oscillations begin when the parameter n exceeds a critical value.")
    print(f"This critical value, n_critical, is found by solving the Hopf bifurcation condition.")
    print(f"The calculated value is:")
    print(f"{n_critical}")

find_critical_n()