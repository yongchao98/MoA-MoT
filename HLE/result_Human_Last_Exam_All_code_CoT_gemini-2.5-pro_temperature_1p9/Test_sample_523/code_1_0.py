import numpy as np
from scipy.optimize import brentq

def find_oscillations_threshold():
    """
    This function calculates the critical Hill coefficient 'n' for which the
    repressilator system starts to exhibit oscillations.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    # We need to find n for which a Hopf bifurcation occurs.
    # The condition for the bifurcation is Re(eigenvalue) = 0, which translates to p - q/2 = 0.
    # where p and q depend on the steady state x_ss, which in turn depends on n.
    
    # 1. Define a function whose root is the steady state x_ss for a given n.
    def steady_state_func(x, n):
        # The equation for the steady state: alpha/(1+x^n) - x - beta/(1+x) = 0
        if x <= 0:
            return float('inf') # Constraint: concentration must be positive
        return alpha / (1 + x**n) - x - beta / (1 + x)

    # 2. Define a function that calculates the bifurcation criterion S = p - q/2.
    #    The root of this function will be the critical value n_c.
    def bifurcation_func(n):
        # For a given n, first find the corresponding steady state x_ss.
        # We need to provide a search interval [a, b] for x_ss where f(a)*f(b) < 0.
        # f(0.1) is positive, f(alpha) is negative.
        try:
            x_ss = brentq(steady_state_func, 1e-3, alpha, args=(n,))
        except ValueError:
            # If no root is found, return a large number to guide the solver.
            return 1e9
        
        # Now, calculate p and q at this steady state.
        p = -1 + beta / (1 + x_ss)**2
        q = -alpha * n * x_ss**(n - 1) / (1 + x_ss**n)**2
        
        # The bifurcation criterion
        S = p - q / 2
        return S

    # 3. Solve for the critical n (n_c) where bifurcation_func(n) = 0.
    #    We need to find a search interval for n. Let's test a few values.
    #    n=2 -> S < 0 (stable)
    #    n=5 -> S > 0 (unstable)
    #    So, the root n_c lies between 2 and 5.
    try:
        n_critical = brentq(bifurcation_func, 2.0, 5.0)
    except ValueError:
        print("Could not find the critical value of n in the specified range.")
        return

    # 4. Find the corresponding steady state, p, and q at this critical point.
    x_ss_critical = brentq(steady_state_func, 1e-3, alpha, args=(n_critical,))
    p_critical = -1 + beta / (1 + x_ss_critical)**2
    q_critical = -alpha * n_critical * x_ss_critical**(n_critical - 1) / (1 + x_ss_critical**n_critical)**2

    # 5. Print the results.
    print("The system exhibits oscillations for values of n greater than a critical value, n_c.")
    print(f"This critical value n_c is found where the system undergoes a Hopf bifurcation.")
    print(f"Calculated critical Hill coefficient, n_c = {n_critical:.4f}\n")
    
    print(f"At this bifurcation point (n = {n_critical:.4f}):")
    print(f"  - The steady state is x_ss = {x_ss_critical:.4f}")
    print(f"  - The Jacobian parameters are p = {p_critical:.4f} and q = {q_critical:.4f}")

    print("\nThe final equation for the bifurcation is p - q/2 = 0.")
    print("Substituting the calculated values:")
    print(f"{p_critical:.4f} - ({q_critical:.4f}) / 2 = {p_critical - q_critical / 2:.4f}")
    
    print(f"\nTherefore, the system exhibits oscillations for n > {n_critical:.4f}")

if __name__ == '__main__':
    find_oscillations_threshold()