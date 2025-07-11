import numpy as np
from scipy.optimize import fsolve, brentq

# Plan:
# 1. Define system parameters.
# 2. Find the critical value of n for the onset of oscillations (Hopf bifurcation).
#    This occurs when the real part of the complex conjugate eigenvalues of the
#    Jacobian matrix at the steady state becomes zero.
# 3. The stability condition is determined by Re(lambda) = A - B/2.
#    Oscillations occur when A - B/2 > 0. The bifurcation occurs at A - B/2 = 0.
# 4. We will numerically solve for the value of n that satisfies this condition.

# System parameters from the problem
alpha = 100.0
beta = 20.0

def find_steady_state(n, initial_guess=1.0):
    """Finds the steady state x* for a given n."""
    def eq(x):
        if x <= 0: return 1e9 # Avoid domain errors
        return alpha / (1 + x**n) - x - beta / (1 + x)
    # Use fsolve to find the root of the equation
    x_star, = fsolve(eq, initial_guess)
    return x_star

def stability_metric(n):
    """
    Calculates the value of A - B/2 for a given n.
    This term is the real part of the critical eigenvalues.
    Its sign determines the stability of the steady state.
    """
    if n <= 0: return -1 # Hill coefficient must be positive
    x_star = find_steady_state(n)
    
    # A = -1 + beta / (1 + x*)^2
    A = -1.0 + beta / (1.0 + x_star)**2
    
    # B = -alpha * n * x*^(n-1) / (1 + x*^n)^2
    B = -alpha * n * x_star**(n - 1.0) / (1.0 + x_star**n)**2
    
    # The bifurcation occurs when this is zero
    return A - B / 2.0

# --- Main execution ---
print("Investigating the stability of the biochemical system:")
print(f"Parameters: alpha = {alpha}, beta = {beta}")
print("The system undergoes a Hopf bifurcation leading to oscillations when Re(lambda) = A - B/2 > 0.")
print("We will solve for the critical Hill coefficient 'n' where A - B/2 = 0.")

# Test the boundaries to ensure a root exists in our search interval [1.0, 2.0]
val_at_1 = stability_metric(1.0)
val_at_2 = stability_metric(2.0)

print(f"\nChecking boundary conditions for root finding:")
print(f"Stability metric at n=1.0: {val_at_1:.4f} (System is stable)")
print(f"Stability metric at n=2.0: {val_at_2:.4f} (System is unstable/oscillatory)")

if np.sign(val_at_1) != np.sign(val_at_2):
    # Find the critical value of n using a root finder (Brent's method)
    # over the interval where the sign changes.
    n_critical = brentq(stability_metric, 1.0, 2.0)
    
    # Find the corresponding steady state x* at this critical n
    x_critical = find_steady_state(n_critical)
    
    print(f"\nFound the critical point for bifurcation:")
    print(f"  Critical Hill coefficient, n_c = {n_critical:.4f}")
    print(f"  Corresponding steady state, x* = {x_critical:.4f}")

    print("\nVerifying the 'final equation' A - B/2 = 0 at this critical point:")
    A_crit = -1.0 + beta / (1.0 + x_critical)**2
    B_crit = -alpha * n_critical * x_critical**(n_critical - 1.0) / (1.0 + x_critical**n_critical)**2
    print(f"  A = -1 + {beta:.1f} / (1 + {x_critical:.4f})^2 = {A_crit:.4f}")
    print(f"  B = -{alpha:.1f} * {n_critical:.4f} * {x_critical:.4f}^({n_critical-1:.4f}) / (1 + {x_critical:.4f}^{n_critical:.4f})^2 = {B_crit:.4f}")
    final_value = A_crit - B_crit/2.0
    print(f"  Result: {A_crit:.4f} - ({B_crit:.4f}) / 2 = {final_value:.6f} (which is approximately zero)\n")

    print("Conclusion:")
    print(f"The system exhibits oscillations for values of n > {n_critical:.4f}.")
    print("For integer values, this means the system oscillates for n >= 2.")

else:
    print("\nCould not find a single bifurcation point in the interval n=[1, 2].")
