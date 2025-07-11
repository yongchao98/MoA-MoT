import numpy as np
from scipy.optimize import fsolve, brentq

def solve_repressilator_bifurcation():
    """
    This script calculates the critical value of the Hill coefficient 'n'
    for the onset of oscillations in a repressilator system.

    The system dynamics are given by:
    dx/dt = alpha / (1 + z^n) - x - beta / (1 + x)
    dy/dt = alpha / (1 + x^n) - y - beta / (1 + y)
    dz/dt = alpha / (1 + y^n) - z - beta / (1 + z)

    Oscillations begin when the symmetric steady state becomes unstable
    via a Hopf bifurcation.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0
    
    print("Analyzing the system with the following parameters:")
    print(f"alpha = {alpha}")
    print(f"beta = {beta}")
    
    def find_steady_state(n, alpha_val=alpha, beta_val=beta):
        """
        Numerically finds the symmetric steady state x* for a given n.
        The steady state equation is: x + beta/(1+x) - alpha/(1+x^n) = 0.
        """
        # Define the function whose root we want to find.
        # We are looking for a positive root x > 0.
        func = lambda x: x[0] + beta_val / (1 + x[0]) - alpha_val / (1 + x[0]**n) if x[0] > 0 else 1e9

        # Provide a reasonable initial guess.
        x_initial_guess = [alpha_val**(1.0 / (n + 1.0))]

        # Use fsolve to find the root.
        x_ss, _, ier, _ = fsolve(func, x_initial_guess, full_output=True)
        
        if ier == 1:
            return x_ss[0]
        else:
            raise RuntimeError(f"Could not find a steady state for n={n}.")

    def bifurcation_condition(n, alpha_val=alpha, beta_val=beta):
        """
        Calculates the real part of the critical eigenvalues.
        The bifurcation occurs when this function is zero.
        Re(lambda) = A - C/2.
        """
        if n <= 0:
            return -1.0 # Stable for non-physical n

        try:
            x_ss = find_steady_state(n, alpha_val, beta_val)

            # Calculate the Jacobian elements at the steady state
            A = -1.0 + beta_val / (1.0 + x_ss)**2
            C = -alpha_val * n * x_ss**(n - 1.0) / (1.0 + x_ss**n)**2

            # The real part of the complex eigenvalues
            real_part = A - C / 2.0
            return real_part
        except RuntimeError:
            return -1.0 # Assume stable if steady state not found

    # Find the critical value of n where the bifurcation occurs.
    # We need to find n such that bifurcation_condition(n) = 0.
    # We bracket the root: a quick test shows it's between n=0.1 (stable) and n=0.5 (unstable).
    try:
        n_critical = brentq(bifurcation_condition, 0.1, 0.5)
    except ValueError:
        print("Error: Could not bracket the root. The function may not cross zero in the given interval.")
        return

    print("\n--- Calculation Result ---")
    print("The final equation for the bifurcation is Re(λ) = A - C/2 = 0, where A and C depend on n.")
    print(f"The critical value of n at which the system starts to oscillate is: {n_critical:.4f}")
    print("For values of n greater than this critical value, the system will exhibit oscillations.")
    
    # Final answer in the required format
    print(f"\n<<<{n_critical}>>>")

# Execute the function
solve_repressilator_bifurcation()