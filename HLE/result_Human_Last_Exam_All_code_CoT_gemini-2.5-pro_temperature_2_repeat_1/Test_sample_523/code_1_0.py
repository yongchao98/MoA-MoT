import numpy as np
from scipy.optimize import fsolve, brentq

def find_oscillation_condition():
    """
    This function analyzes the stability of the given biochemical system
    to find the values of n for which it exhibits oscillations.
    """
    # System parameters
    alpha = 100
    beta = 20

    print("To find the conditions for oscillation, we analyze the stability of the system's symmetric fixed point (x*, x*, x*).")
    print("Oscillations arise from a Hopf bifurcation when the fixed point loses stability.")
    print("\nFirst, the fixed point x* must satisfy the equation:")
    print("α / (1 + x*^n) - x* - β / (1 + x*) = 0")
    print(f"With α={alpha} and β={beta}, this is:")
    print(f"{alpha} / (1 + x*^n) - x* - {beta} / (1 + x*) = 0")

    print("\nSecond, the stability is determined by the eigenvalues of the system's Jacobian matrix at the fixed point.")
    print("A Hopf bifurcation occurs when the real part of a pair of complex conjugate eigenvalues becomes zero.")
    print("This leads to the following stability condition, which must be zero at the bifurcation point:")
    print("-1 + β/((1+x*)²) + (α*n*x*^(n-1))/(2*(1+x*^n)²) = 0")
    print("Plugging in the parameter values, we need to solve for the critical n (n_c):")
    print(f"-1 + {beta}/((1+x*)²) + ({alpha}*n_c*x*^(n_c-1))/(2*(1+x*^(n_c))²) = 0")
    
    # --- Numerical Solution ---
    
    # Function for the fixed point equation: f(x, n) = 0
    def fixed_point_eq(x, n):
        # We need to handle potential negative values of x during root finding
        if x < 0:
            return np.inf
        return alpha / (1 + x**n) - x - beta / (1 + x)

    # Function that returns the stability condition value.
    # We want to find the root of this function of n.
    def stability_bifurcation_eq(n):
        # For a given n, first find the corresponding fixed point x*
        # The equation for x* is monotonically decreasing, so fsolve is reliable
        # with a reasonable initial guess.
        initial_guess_x = 3.0  # Based on preliminary analysis
        x_star, = fsolve(fixed_point_eq, initial_guess_x, args=(n,))
        
        # Now, calculate the value of the stability condition
        # It's positive for unstable (oscillating) systems
        return -1 + beta / ((1 + x_star)**2) + (alpha * n * x_star**(n - 1)) / (2 * (1 + x_star**n)**2)

    # We need to find the root n_c of stability_bifurcation_eq(n) = 0.
    # We will use a bracketing root-finding algorithm (brentq), so we need
    # to find an interval [n_low, n_high] where the function changes sign.
    n_low = 1.0
    n_high = 3.0 # Let's test n=1 and n=3 to find a bracket
    
    val_low = stability_bifurcation_eq(n_low)
    val_high = stability_bifurcation_eq(n_high)

    print(f"\nNumerically evaluating the stability condition:")
    print(f"For n = {n_low}, the stability value is {val_low:.4f} (stable).")
    print(f"For n = {n_high}, the stability value is {val_high:.4f} (unstable).")
    
    # Since the signs are different, a root exists in the interval [1.0, 3.0].
    # We can now find the critical value n_c precisely.
    if np.sign(val_low) != np.sign(val_high):
        n_critical = brentq(stability_bifurcation_eq, n_low, n_high)
        print("\nSolving for the critical value n_c where the stability condition equals zero, we get:")
        print(f"n_c = {n_critical:.4f}")
        print("\nTherefore, the system exhibits oscillations for values of n > n_c.")
        print(f"The condition for oscillations is n > {n_critical:.4f}.")
        print(f"The smallest integer value for n that produces oscillations is {int(np.ceil(n_critical))}.")
    else:
        print("\nCould not find a bracket for the root. The stability does not change in the tested range.")

# Execute the analysis
find_oscillation_condition()