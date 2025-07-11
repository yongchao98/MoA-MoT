import numpy as np
from scipy.optimize import fsolve

def solve_oscillation_condition():
    """
    This function encapsulates the full analysis to find the critical value of 'n'
    for the onset of oscillations in the given biochemical system.
    """
    
    # --- System Parameters ---
    alpha = 100.0
    beta = 20.0

    print("Analyzing the biochemical system:")
    print("dx/dt = alpha / (1 + z^n) - x - beta / (1 + x)")
    print("dy/dt = alpha / (1 + x^n) - y - beta / (1 + y)")
    print("dz/dt = alpha / (1 + y^n) - z - beta / (1 + z)")
    print(f"with alpha = {alpha} and beta = {beta}\n")

    # --- Theoretical Framework ---
    # The onset of oscillations corresponds to a Hopf bifurcation.
    # This occurs when the real part of a pair of complex conjugate eigenvalues
    # of the system's Jacobian matrix at the fixed point becomes zero.
    # For this system, the condition simplifies to: A - B/2 = 0
    # where A and B are terms derived from the Jacobian evaluated at the symmetric
    # fixed point (x*, x*, x*).
    
    def bifurcation_condition(n, alpha_val, beta_val):
        """
        Calculates the stability metric (A - B/2) for a given n.
        The root of this function is the critical n for the Hopf bifurcation.
        """
        if n <= 0: return 1e6

        fixed_point_eq = lambda x: x + beta_val / (1.0 + x) - alpha_val / (1.0 + np.power(x, n))
        
        try:
            # Solve for the symmetric fixed point x*
            x_star, = fsolve(fixed_point_eq, 1.0, xtol=1e-9, maxfev=500)
            if x_star < 0: return 1e6
        except Exception:
            return 1e6 # Return a large value if solver fails

        # Calculate Jacobian-derived terms A and B
        A = -1.0 + beta_val / (1.0 + x_star)**2
        
        try:
            x_pow_n = np.power(x_star, n)
            x_pow_n_minus_1 = np.power(x_star, n - 1)
        except OverflowError:
            return 1e6
            
        B_numerator = -alpha_val * n * x_pow_n_minus_1
        B_denominator = (1.0 + x_pow_n)**2
        
        B = B_numerator / B_denominator if B_denominator != 0 else 0.0

        return A - B / 2.0

    # --- Numerical Solution ---
    print("Searching for the critical value of 'n' where oscillations begin...")
    
    # Find the root of the bifurcation condition to get n_critical
    try:
        n_critical, = fsolve(bifurcation_condition, x0=1.5, args=(alpha, beta), xtol=1e-9)
    except Exception as e:
        print(f"Failed to find the solution. Error: {e}")
        return None

    # --- Results ---
    print(f"\nCalculation complete.")
    print(f"The critical value for the onset of oscillations is n_critical = {n_critical:.4f}")

    # Re-calculate values at n_critical to display them as requested
    x_star_crit, = fsolve(lambda x: x + beta / (1.0 + x) - alpha / (1.0 + np.power(x, n_critical)), 1.0)
    A_crit = -1.0 + beta / (1.0 + x_star_crit)**2
    x_pow_n_crit = np.power(x_star_crit, n_critical)
    x_pow_n_minus_1_crit = np.power(x_star_crit, n_critical - 1)
    B_crit = (-alpha * n_critical * x_pow_n_minus_1_crit) / (1.0 + x_pow_n_crit)**2

    print("\nAt the point of bifurcation (n = n_critical):")
    print(f"  The fixed point of the system is x* = y* = z* = {x_star_crit:.4f}")
    
    print(f"\n  The equation for the bifurcation is: A - B/2 = 0")
    print("  The values of the components of this equation are:")
    print(f"    A = -1 + {beta} / (1 + {x_star_crit:.4f})^2 = {A_crit:.4f}")
    print(f"    B = (-{alpha} * {n_critical:.4f} * {x_star_crit:.4f}^({n_critical:.4f}-1)) / (1 + {x_star_crit:.4f}^{n_critical:.4f})^2 = {B_crit:.4f}")
    
    print(f"\n  Substituting these values into the equation:")
    print(f"  {A_crit:.4f} - ({B_crit:.4f} / 2) = {A_crit - B_crit / 2.0:.6f} (which is effectively zero, as expected)")
    
    print("\nFor any value of n greater than this critical value, the system will exhibit oscillations.")
    print(f"Conclusion: Oscillations occur for n > {n_critical:.4f}")
    
    return n_critical

if __name__ == '__main__':
    solve_oscillation_condition()