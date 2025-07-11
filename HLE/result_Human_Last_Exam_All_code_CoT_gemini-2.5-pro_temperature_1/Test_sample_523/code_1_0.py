import numpy as np
from scipy.optimize import fsolve, brentq

def solve_system():
    """
    Analyzes the biochemical system to find the values of n for which it oscillates.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    def steady_state_equation(x, n):
        """
        The equation for the symmetric steady state (x=y=z).
        We need to find the root of this function, which is x*.
        f(x*) = 0
        """
        # To avoid overflow with large x or n, we check the magnitude.
        # However, for this problem, x* is expected to be in a reasonable range.
        if x <= 0:
            return 1e12 # Return a large number for non-physical values
        return alpha / (1.0 + x**n) - x - beta / (1.0 + x)

    def bifurcation_condition(n):
        """
        The Hopf bifurcation condition, A - B/2 = 0.
        Oscillations occur when this function is > 0.
        We find the root of this function to get the critical n.
        """
        # First, find the steady state x* for the given n.
        # We use an initial guess of 3.0, as our initial analysis showed the root
        # is likely in that region.
        # The args parameter passes the additional arguments (n) to steady_state_equation.
        x_star_result = fsolve(steady_state_equation, 3.0, args=(n,))
        x_star = x_star_result[0]
        
        # If the solver fails or finds a non-physical root, return a value
        # that won't be a root.
        if x_star <= 0:
            return -1.0

        # Calculate the terms of the Jacobian matrix, A and B
        # A = -1 + beta / (1 + x*)^2
        # B = -alpha * n * x*^(n-1) / (1 + x*^n)^2
        A = -1.0 + beta / (1.0 + x_star)**2
        B = -alpha * n * x_star**(n - 1.0) / (1.0 + x_star**n)**2
        
        # The condition for the onset of oscillations (Hopf bifurcation) is Re(lambda) = 0.
        # This corresponds to A - B/2 = 0.
        return A - B/2

    # We need to find the root of the bifurcation_condition function.
    # From a preliminary check, we know the root lies between n=1 and n=2.
    # We use the Brent method for secure root finding in a given bracket.
    try:
        n_critical = brentq(bifurcation_condition, 1.0, 2.0)
    except ValueError:
        print("Could not find a root in the specified interval [1.0, 2.0].")
        print("The bifurcation condition may not change sign in this interval.")
        return

    # Now we have the critical value of n. Let's present the results.
    print(f"The system is given by the equations:")
    print(f"dx/dt = {alpha}/(1+z^n) - x - {beta}/(1+x)")
    print(f"dy/dt = {alpha}/(1+x^n) - y - {beta}/(1+y)")
    print(f"dz/dt = {alpha}/(1+y^n) - z - {beta}/(1+z)\n")
    
    print(f"Oscillations occur when the system's symmetric steady state becomes unstable.")
    print(f"This happens via a Hopf bifurcation when the cooperativity coefficient 'n' exceeds a critical value.")
    print(f"\nThe critical value of n is found to be: {n_critical:.4f}")
    print(f"Therefore, the system exhibits oscillations for n > {n_critical:.4f}.\n")

    # As requested, output the numbers in the final equation.
    # The "final equation" is the bifurcation condition A - B/2 = 0.
    # Let's calculate the values at the critical point.
    x_star_critical = fsolve(steady_state_equation, 3.0, args=(n_critical,))[0]
    A_critical = -1.0 + beta / (1.0 + x_star_critical)**2
    B_critical = -alpha * n_critical * x_star_critical**(n_critical - 1.0) / (1.0 + x_star_critical**n_critical)**2
    
    print("This critical value is determined by the stability condition A - B/2 = 0.")
    print("At the point of bifurcation:")
    print(f"  n_critical = {n_critical:.4f}")
    print(f"  Steady state x* = {x_star_critical:.4f}")
    print(f"  Jacobian term A = -1 + {beta} / (1 + {x_star_critical:.4f})^2 = {A_critical:.4f}")
    print(f"  Jacobian term B = -{alpha} * {n_critical:.4f} * {x_star_critical:.4f}^({n_critical:.4f}-1) / (1 + {x_star_critical:.4f}^{n_critical:.4f})^2 = {B_critical:.4f}")
    print(f"\nPlugging these into the equation A - B/2:")
    print(f"  {A_critical:.4f} - ({B_critical:.4f}) / 2 = {A_critical - B_critical/2:.6f} (which is approximately 0 as expected)")

solve_system()