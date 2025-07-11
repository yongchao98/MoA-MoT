import numpy as np
from scipy.optimize import fsolve, brentq

def solve_for_oscillations():
    """
    This script analyzes the provided system of differential equations to find the
    parameter values of 'n' that lead to oscillations.
    """
    # Define the system parameters from the problem description
    alpha = 100.0
    beta = 20.0

    # Define the equation for the fixed point x*.
    # The fixed point is a root of f(x) = x + beta/(1+x) - alpha/(1+x^n) = 0.
    def fixed_point_eq(x, n, alpha, beta):
        # We only consider positive concentrations x > 0.
        if x <= 0:
            # Return a large number to guide the solver away from non-physical solutions.
            return np.inf
        return x + beta / (1.0 + x) - alpha / (1.0 + x**n)

    # Define the function for the Hopf bifurcation condition.
    # Oscillations begin when this function's value crosses from negative to positive.
    # The condition is derived from the eigenvalues of the system's Jacobian matrix.
    def bifurcation_condition(n, alpha, beta):
        """
        Calculates the value of the Hopf bifurcation condition, Re(lambda) = A - B/2.
        The root of this function gives the critical value of n.
        """
        # First, find the fixed point x* for the given n.
        # An initial guess of x0=3.0 is robust for the expected range of n.
        try:
            x_star = fsolve(fixed_point_eq, x0=3.0, args=(n, alpha, beta))[0]
            if x_star <= 0:
                # If a non-physical fixed point is found, we are not at the bifurcation.
                return -1.0
        except:
            # If the solver fails, assume we are not at the bifurcation.
            return -1.0

        # Calculate the terms A and B from the Jacobian matrix at the fixed point.
        # A = -1 + beta / (1 + x*)^2
        # B = -alpha * n * x*^(n-1) / (1 + x*^n)^2
        A = -1.0 + beta / (1.0 + x_star)**2
        B = -alpha * n * x_star**(n - 1) / (1.0 + x_star**n)**2
        
        # The real part of the complex eigenvalues is (A - B/2).
        # The Hopf bifurcation occurs when this term is zero.
        return A - B / 2.0

    try:
        # We search for the root of the bifurcation condition in a reasonable interval for n.
        # Manual checks show the root lies between n=1 and n=2.
        # The brentq algorithm is efficient for finding a root within a bracket.
        n_critical = brentq(bifurcation_condition, a=1.0, b=2.0, args=(alpha, beta))
        
        # For this critical n, we can find the corresponding fixed point concentration x*
        x_star_critical = fsolve(fixed_point_eq, x0=3.0, args=(n_critical, alpha, beta))[0]

        print("For the system with parameters:")
        print(f"alpha = {int(alpha)}")
        print(f"beta = {int(beta)}")
        
        print("\nOscillations emerge via a Hopf bifurcation.")
        print("The critical value of the parameter 'n' at which this occurs is found by solving the bifurcation condition.")
        
        print(f"\nFinal Equation (Bifurcation Condition): A - B/2 = 0")
        print(f"Calculated critical value of n: {n_critical:.4f}")
        print(f"Corresponding steady-state concentration x*: {x_star_critical:.4f}")

        print(f"\nConclusion: The system exhibits oscillations for values of n > {n_critical:.4f}.")

        # The final answer format requires a single value. We return the critical n.
        return n_critical

    except ValueError:
        print("Error: The critical value for n could not be found in the specified interval [1, 2].")
        return None

# Execute the analysis and store the result
critical_n_value = solve_for_oscillations()

# The final answer is the critical value of n
# print(f"\n<<< {critical_n_value:.4f} >>>")