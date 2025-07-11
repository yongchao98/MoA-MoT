import numpy as np
from scipy.optimize import root_scalar

def solve_for_n():
    """
    This script finds the critical value of the parameter 'n' for which the
    given biochemical system starts to exhibit oscillations.
    """
    # Parameters from the problem statement
    alpha = 100
    beta = 20

    def steady_state_func(x, n):
        """
        Defines the steady-state equation f(x, n) = 0.
        x_s + beta / (1 + x_s) - alpha / (1 + x_s^n) = 0
        """
        if x <= 0:
            return -np.inf
        try:
            return x + beta / (1 + x) - alpha / (1 + x**n)
        except (OverflowError, ValueError):
            return -np.inf

    def bifurcation_condition_func(n):
        """
        Calculates the value of Re(lambda_2,3) = a - b/2 for a given n.
        A root of this function corresponds to the Hopf bifurcation point.
        """
        # First, find the steady state x_s corresponding to this n
        try:
            sol = root_scalar(steady_state_func, args=(n,), bracket=[0.1, alpha], method='brentq')
            x_s = sol.root
        except ValueError:
            # If finding the steady state fails, this 'n' is not valid
            return None

        # Now, calculate a - b/2 at (x_s, n)
        # a = -1 + beta / (1 + x_s)^2
        term_a = -1 + beta / (1 + x_s)**2

        # b = -alpha * n * x_s^(n-1) / (1 + x_s^n)^2
        try:
            x_s_n = x_s**n
            x_s_n_minus_1 = x_s**(n-1)
            term_b = -alpha * n * x_s_n_minus_1 / (1 + x_s_n)**2
        except (OverflowError, ValueError):
             # For large n, x_s**n can overflow. Handle gracefully.
            if n * np.log(x_s) > 50: # log(very_large_number)
                term_b = 0.0 # Denominator becomes infinite
            else: # Other numerical issues
                return None


        return term_a - term_b / 2

    # --- Main Calculation ---
    # Find the root of bifurcation_condition_func to get the critical n.
    # From a preliminary analysis, we know the root is between n=1 and n=2.
    try:
        critical_n_sol = root_scalar(bifurcation_condition_func, bracket=[1.0, 2.0], method='brentq')
        if critical_n_sol.converged:
            n_c = critical_n_sol.root
            print("To find the values of 'n' for which the system oscillates, we analyze the system's stability.")
            print("Oscillations begin at a Hopf bifurcation, where a steady state becomes unstable.")
            print(f"\nThe system parameters are:")
            print(f"  alpha = {alpha}")
            print(f"  beta = {beta}")
            print("\nThe bifurcation occurs at a critical value of 'n' (n_c).")
            print(f"Numerically solving for this value, we find:")
            print(f"  n_c ≈ {n_c:.4f}")

            # For oscillations, Re(lambda) must be positive. This occurs for n > n_c.
            print("\nThe system exhibits oscillations for values of n greater than this critical value.")
            print(f"Therefore, the condition for oscillations is: n > {n_c:.4f}")

        else:
            print("Numerical solver for n_c did not converge.")

    except Exception as e:
        print(f"An error occurred: {e}")

solve_for_n()

# The question asks for the values of n, which is the range where oscillations occur.
# The numeric value that defines this range is the critical n.
<<<1.2556>>>