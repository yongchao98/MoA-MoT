import numpy as np
from scipy.optimize import root, brentq

def solve_repressilator_stability():
    """
    This script finds the critical value of the parameter 'n' for which the
    repressilator system with titration exhibits oscillations.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    # --- Step 1: Define function to find the fixed point ---

    def fixed_point_equation(x, n, alpha, beta):
        """Equation for the symmetric fixed point: f(x, n) = 0."""
        x_val = x[0]
        if x_val <= 0:
            return [1e9]  # Return a large number for invalid domain
        return [alpha / (1 + x_val**n) - x_val - beta / (1 + x_val)]

    def find_fixed_point(n, alpha, beta):
        """Finds the fixed point x* for a given n using a robust solver."""
        # Initial guess for the root finder.
        # The fixed point value changes with n, so we try a few robust guesses if one fails.
        initial_guesses = [4.0, 2.5, 8.0]
        for guess in initial_guesses:
            solution = root(fixed_point_equation, [guess], args=(n, alpha, beta), method='hybr')
            if solution.success:
                return solution.x[0]
        raise RuntimeError(f"Failed to find a fixed point for n={n}")

    # --- Step 2: Define the Hopf bifurcation condition ---

    def hopf_bifurcation_condition(n, alpha, beta):
        """
        Calculates the value of the Hopf bifurcation criterion, A - B/2.
        The system oscillates when this value is > 0. We search for its root.
        """
        if n <= 0:
            return 1e9  # n must be positive

        try:
            x_star = find_fixed_point(n, alpha, beta)
            
            # A = -1 + beta / (1 + x*)^2
            A = -1.0 + beta / (1.0 + x_star)**2
            
            # B = -alpha * n * x*^(n-1) / (1 + x*^n)^2
            B = -alpha * n * x_star**(n - 1) / (1.0 + x_star**n)**2
            
            # The condition for the onset of oscillations is Re(lambda) = A - B/2 = 0
            return A - B / 2.0
        except RuntimeError:
            # If fixed point finding fails, return a large number
            return 1e9

    # --- Step 3: Find the critical value of n ---

    # Use a root finder to find n_critical where the bifurcation condition is zero.
    # We search in a reasonable interval, e.g., [1.0, 3.0].
    try:
        n_critical = brentq(hopf_bifurcation_condition, a=1.0, b=3.0, args=(alpha, beta))
    except ValueError:
        print("Could not find the critical value of n in the specified interval.")
        return

    # --- Step 4: Output the results ---
    
    print(f"The system is described by the equations:")
    print(f"dx/dt = alpha / (1 + z^n) - x - beta / (1 + x)")
    print(f"dy/dt = alpha / (1 + x^n) - y - beta / (1 + y)")
    print(f"dz/dt = alpha / (1 + y^n) - z - beta / (1 + z)")
    print("\nWith parameter values:")
    print(f"alpha = {alpha}")
    print(f"beta = {beta}")
    
    print("\nOscillations occur when the system's fixed point becomes unstable via a Hopf bifurcation.")
    print("This happens for values of n greater than a critical value, n_critical.")
    
    print(f"\nThe calculated critical value is:")
    print(f"n_critical = {n_critical}")
    
    # For verification, let's show the values at the critical point
    x_star_critical = find_fixed_point(n_critical, alpha, beta)
    A_critical = -1.0 + beta / (1.0 + x_star_critical)**2
    B_critical = -alpha * n_critical * x_star_critical**(n_critical - 1) / (1.0 + x_star_critical**n_critical)**2
    
    print("\nAt this critical point, the final equation for the bifurcation 'A - B/2 = 0' is satisfied:")
    print(f"Fixed point x* = {x_star_critical}")
    print(f"A = {A_critical}")
    print(f"B = {B_critical}")
    print(f"A - B/2 = {A_critical - B_critical / 2.0:.2e}")
    
    print("\nTherefore, the system exhibits oscillations for n >", n_critical)
    
    # Return the final numerical answer for the <<<>>> tag
    return n_critical

# Execute the analysis and capture the final answer
final_answer = solve_repressilator_stability()
print(f"\n<<<{final_answer}>>>")
