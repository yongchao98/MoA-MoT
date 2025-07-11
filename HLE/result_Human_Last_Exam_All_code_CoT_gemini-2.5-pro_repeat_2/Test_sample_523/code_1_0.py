import numpy as np
from scipy.optimize import fsolve

def solve_for_oscillations():
    """
    This function determines the range of the parameter 'n' for which the
    given biochemical system exhibits oscillations.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    # The system of equations is:
    # dx/dt = alpha / (1 + z^n) - x - beta / (1 + x)
    # dy/dt = alpha / (1 + x^n) - y - beta / (1 + y)
    # dz/dt = alpha / (1 + y^n) - z - beta / (1 + z)

    # --- Step 1: Define function to find the steady state x* ---
    # At steady state x=y=z=x*, the equation becomes:
    # alpha / (1 + x*^n) - x* - beta / (1 + x*) = 0
    def steady_state_eq(x, n, alpha_val, beta_val):
        """The steady-state equation we need to find the root of."""
        return alpha_val / (1 + x**n) - x - beta_val / (1 + x)

    def find_steady_state(n, alpha_val, beta_val):
        """Numerically solves for the steady state x* for a given n."""
        # We use an initial guess of 1.0 for x*
        x_star, = fsolve(steady_state_eq, 1.0, args=(n, alpha_val, beta_val))
        return x_star

    # --- Step 2: Define the Hopf bifurcation condition ---
    # A Hopf bifurcation occurs when the real part of the complex eigenvalues of
    # the Jacobian matrix is zero. For this system, this condition is:
    # A - B/2 = 0
    # where:
    # A = -1 + beta / (1 + x*)^2
    # B = -alpha * n * x*^(n-1) / (1 + x*^n)^2
    def hopf_condition(n, alpha_val, beta_val):
        """
        Represents the Hopf condition (Re(lambda_complex) = 0). We will find
        the value of 'n' that makes this function equal to zero.
        """
        if n <= 0: # Hill coefficient n must be positive
            return np.inf

        # For a given n, first find the corresponding steady state x*
        x_star = find_steady_state(n, alpha_val, beta_val)
        
        if x_star <= 0: # Concentration must be positive
            return np.inf

        A = -1.0 + beta_val / (1.0 + x_star)**2
        B = -alpha_val * n * x_star**(n - 1) / (1.0 + x_star**n)**2
        
        return A - B / 2.0

    # --- Step 3: Solve for the critical value of n ---
    try:
        # We use fsolve to find the root of the hopf_condition function,
        # which gives the critical value of n. We start the search from n=2.0.
        n_critical = fsolve(hopf_condition, 2.0, args=(alpha, beta))[0]
        
        # --- Step 4: Output the results ---
        print("The system of equations is given by:")
        print(f"dx/dt = {alpha}/{'(1+z^n)'} - x - {beta}/{'(1+x)'}")
        print(f"dy/dt = {alpha}/{'(1+x^n)'} - y - {beta}/{'(1+y)'}")
        print(f"dz/dt = {alpha}/{'(1+y^n)'} - z - {beta}/{'(1+z)'}")
        print("\nFor the given parameters:")
        print(f"alpha = {int(alpha)}")
        print(f"beta = {int(beta)}")
        
        print("\nOscillations occur when the system's steady state becomes unstable via a Hopf bifurcation.")
        print(f"This bifurcation happens at a critical value of the Hill coefficient, n.")
        print(f"The calculated critical value is n_c = {n_critical:.4f}.")
        print("\nThe system exhibits oscillations for values of n greater than this critical value.")
        print(f"Therefore, the condition for oscillation is: n > {n_critical:.4f}")

    except Exception as e:
        print(f"An error occurred during the calculation: {e}")

if __name__ == "__main__":
    solve_for_oscillations()