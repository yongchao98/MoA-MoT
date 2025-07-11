import numpy as np
from scipy.optimize import fsolve

def solve_repressilator_bifurcation():
    """
    This script determines the critical value of the Hill coefficient 'n'
    for which a symmetric three-gene repressilator system starts to oscillate.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    def find_steady_state(n_val, alpha_val, beta_val):
        """
        Numerically finds the symmetric steady-state concentration 's' for a given 'n'.
        We need to solve the equation: alpha / (1 + s^n) - s - beta / (1 + s) = 0.
        """
        # The function whose root we want to find.
        # s is passed as a single-element array by fsolve.
        func_to_solve = lambda s: alpha_val / (1 + s[0]**n_val) - s[0] - beta_val / (1 + s[0])
        
        # An initial guess for s. Based on analysis, the solution is expected
        # to be near the minimum of g(s) = s + beta/(1+s), which occurs at s=sqrt(beta)-1.
        initial_guess_s = [np.sqrt(beta_val) - 1]
        s_solution, = fsolve(func_to_solve, initial_guess_s)
        return s_solution

    def hopf_bifurcation_condition(n_val, alpha_val, beta_val):
        """
        Calculates the value of the Hopf bifurcation condition, A - B/2.
        We are looking for the value of n that makes this function zero.
        n_val is passed as a single-element array by fsolve.
        """
        n = n_val[0]
        s = find_steady_state(n, alpha_val, beta_val)

        # A concentration must be positive for a physically realistic solution.
        if s <= 0:
            return np.inf

        # A and B are terms from the Jacobian matrix evaluated at the steady state (s,s,s).
        # J = [[A, 0, B], [B, A, 0], [0, B, A]]
        A = -1.0 + beta_val / (1.0 + s)**2
        B = -alpha_val * n * s**(n - 1.0) / (1.0 + s**n)**2

        # The bifurcation occurs when the real part of the complex eigenvalues (A - B/2) is zero.
        return A - B / 2.0

    # Find the critical value of n by finding the root of the Hopf condition function.
    # From preliminary analysis, the root is expected to be between 1 and 2.
    initial_guess_n = [1.5]
    n_critical, = fsolve(hopf_bifurcation_condition, initial_guess_n, args=(alpha, beta))

    # --- Output the Results ---
    print(f"The system exhibits oscillations for values of the Hill coefficient 'n' greater than a critical value.")
    print(f"The critical value is n_critical = {n_critical:.4f}\n")
    
    print("This critical value is found at the Hopf bifurcation, where the condition A - B/2 = 0 is met.")
    
    # Calculate and print the values at the critical point for verification
    s_critical = find_steady_state(n_critical, alpha, beta)
    A_critical = -1.0 + beta / (1.0 + s_critical)**2
    B_critical = -alpha * n_critical * s_critical**(n_critical - 1.0) / (1.0 + s_critical**n_critical)**2
    
    print("--- Verification at the Critical Point ---")
    print(f"The final equation to solve for the bifurcation is: A - B/2 = 0")
    print(f"At the calculated critical n = {n_critical:.4f}, the corresponding values are:")
    print(f"  Steady State (s) = {s_critical:.4f}")
    print(f"  Jacobian Term (A) = {A_critical:.4f}")
    print(f"  Jacobian Term (B) = {B_critical:.4f}")
    
    # Show the numbers in the final equation
    print("\nPlugging these numbers into the equation:")
    print(f"  {A_critical:.4f} - ({B_critical:.4f}) / 2 = {A_critical - B_critical / 2.0:.6f}")
    print("This result is approximately zero, confirming our calculation.")

solve_repressilator_bifurcation()
<<<1.3366>>>