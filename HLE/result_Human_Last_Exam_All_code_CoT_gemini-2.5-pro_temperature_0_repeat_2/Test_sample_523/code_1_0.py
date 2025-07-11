import numpy as np
from scipy.optimize import root_scalar

def solve_for_oscillations():
    """
    This function calculates the critical value of 'n' for the onset of oscillations
    in the given biochemical system.
    """
    # Define the given parameters
    alpha = 100.0
    beta = 20.0

    def steady_state_eq(s, n):
        """
        Defines the steady-state equation: s + beta/(1+s) - alpha/(1+s^n) = 0.
        This function is used by a root finder to solve for 's'.
        """
        if s <= 0:
            # Return a large number for invalid 's' to guide the solver
            return 1e10
        return s + beta / (1 + s) - alpha / (1 + s**n)

    def find_s_for_n(n):
        """
        Numerically solves the steady-state equation for 's' for a given 'n'.
        """
        # Use a root-finding algorithm to find 's'. The bracket [0.01, 100] is a safe
        # range where a solution is expected to exist.
        sol = root_scalar(steady_state_eq, args=(n,), bracket=[0.01, 100], method='brentq')
        return sol.root

    def hopf_condition(n):
        """
        Represents the Hopf bifurcation condition, A - B/2 = 0.
        This function is used by a root finder to solve for the critical 'n'.
        """
        # First, find the corresponding steady state 's' for the given 'n'
        s = find_s_for_n(n)

        # Calculate the terms A and B from the Jacobian matrix
        A = -1.0 + beta / (1.0 + s)**2
        B = -alpha * n * s**(n - 1) / (1.0 + s**n)**2

        # The Hopf bifurcation occurs when the real part of the complex eigenvalues is zero
        return A - B / 2.0

    # Find the critical value of 'n' by finding the root of the Hopf condition.
    # From preliminary analysis, we know the root is between 1.1 and 1.2.
    try:
        n_critical_sol = root_scalar(hopf_condition, bracket=[1.1, 1.2], method='brentq')
        n_c = n_critical_sol.root
        s_c = find_s_for_n(n_c)

        print(f"The system exhibits oscillations for values of n greater than a critical value, n_c.")
        print(f"The calculated critical value is n_c = {n_c:.4f}\n")
        
        print("This critical point is where the Hopf bifurcation occurs. We can verify the condition A - B/2 = 0.")
        print("The final equation for the bifurcation is:")
        print("(-1 + beta / (1 + s)^2) - (-alpha * n * s^(n-1) / (1 + s^n)^2) / 2 = 0\n")
        
        print("Plugging in the calculated values at the bifurcation point:")
        A_val = -1.0 + beta / (1.0 + s_c)**2
        B_val = -alpha * n_c * s_c**(n_c - 1) / (1.0 + s_c**n_c)**2
        
        print(f"n = {n_c:.4f}")
        print(f"s = {s_c:.4f}")
        print(f"alpha = {alpha}")
        print(f"beta = {beta}\n")

        print("The equation becomes:")
        print(f"(-1 + {beta} / (1 + {s_c:.4f})^2) - (-{alpha} * {n_c:.4f} * {s_c:.4f}^({n_c:.4f - 1:.4f}) / (1 + {s_c:.4f}^{n_c:.4f})^2) / 2 = 0")
        print(f"({A_val:.4f}) - ({B_val:.4f}) / 2 = {A_val - B_val/2:.2e}")
        
        # The final answer format requires a single value.
        print(f"\nThus, the system exhibits oscillations for n > {n_c:.4f}")
        print(f"<<<{n_c:.4f}>>>")

    except ValueError:
        print("Error: Could not find the critical value of n in the specified bracket.")

# Run the solver
solve_for_oscillations()