import numpy as np
from scipy.optimize import fsolve, brentq

def solve_biochemical_oscillations():
    """
    This script finds the critical value of the parameter 'n' for which the given
    biochemical system starts to exhibit oscillations.
    """
    
    # System parameters
    alpha = 100.0
    beta = 20.0

    def steady_state_equation(x, n, alpha, beta):
        """
        The equation that defines the symmetric steady state x_ss.
        We need to find the root of this function.
        """
        # Concentrations must be non-negative.
        if x < 0:
            return 1e9
        return alpha / (1 + x**n) - x - beta / (1 + x)

    def find_bifurcation_value(n, alpha, beta):
        """
        Calculates the value of the Hopf bifurcation condition, A - B/2.
        The root of this function is the critical value of n.
        """
        # Find the steady state x_ss for the given n using a numerical solver.
        # An initial guess of 3.0 is used based on preliminary analysis.
        try:
            x_ss = fsolve(steady_state_equation, 3.0, args=(n, alpha, beta))[0]
            if x_ss < 0:  # Not a physical solution
                return 1e9  # Return a large number to steer solver away
        except:
            return 1e9 # fsolve might fail for some n

        # Calculate the Jacobian matrix elements A and B at the steady state
        A = -1.0 + beta / ((1 + x_ss)**2)
        B = -alpha * n * x_ss**(n - 1) / ((1 + x_ss**n)**2)

        # The Hopf bifurcation occurs when A - B/2 = 0
        return A - B / 2.0

    # Use a root-finding algorithm (Brent's method) to find the critical n.
    # From manual checks, we know n=1 is stable and n=2 is oscillatory,
    # so the root must lie in the interval [1.0, 2.0].
    try:
        n_critical = brentq(find_bifurcation_value, 1.0, 2.0, args=(alpha, beta))

        # The system exhibits oscillations for n > n_critical.
        print(f"The system exhibits oscillations for values of n > {n_critical:.4f}")
        print("\nThis critical value is where a Hopf bifurcation occurs, leading to oscillatory behavior.")
        
        # To fulfill the requirement of showing the numbers in the final equation,
        # we calculate the values at the bifurcation point.
        x_ss_crit = fsolve(steady_state_equation, 3.0, args=(n_critical, alpha, beta))[0]
        A_crit = -1.0 + beta / ((1 + x_ss_crit)**2)
        B_crit = -alpha * n_critical * x_ss_crit**(n_critical - 1) / ((1 + x_ss_crit**n_critical)**2)
        condition_value = A_crit - B_crit / 2.0
        
        print("\nAt the critical point:")
        print(f"n = {n_critical:.4f}")
        print(f"Steady state x_ss = {x_ss_crit:.4f}")
        print(f"Jacobian element A = {A_crit:.4f}")
        print(f"Jacobian element B = {B_crit:.4f}")
        
        print("\nThe bifurcation condition is A - B/2 = 0. Substituting the values at the critical point:")
        print(f"({A_crit:.4f}) - ({B_crit:.4f}) / 2 = {condition_value:.4f}")

    except ValueError:
        print("Could not find the critical value of n in the specified interval [1.0, 2.0].")

solve_biochemical_oscillations()