import numpy as np
from scipy.optimize import root

def solve_for_n():
    """
    This function finds the critical Hill coefficient 'n' for the onset of oscillations
    in the given biochemical system.
    """
    
    # System parameters as given in the problem
    alpha = 100.0
    beta = 20.0
    
    # We define the system of two nonlinear equations to find the critical point (s_c, n_c).
    # The input `x` is a list or array where x[0] corresponds to s and x[1] corresponds to n.
    def find_bifurcation_point(x):
        s, n = x[0], x[1]
    
        # Return a large error for invalid physical values (s and n must be positive)
        if s <= 0 or n <= 0:
            return [1e9, 1e9]
    
        try:
            # Pre-calculate s^n as it's used multiple times, preventing re-computation
            s_pow_n = s**n
            
            # Equation 1: The fixed point condition
            # s + beta/(1+s) - alpha/(1+s^n) = 0
            eq1 = s + beta / (1 + s) - alpha / (1 + s_pow_n)
    
            # Equation 2: The Hopf bifurcation condition
            # 1 - beta/((1+s)^2) - (alpha * n * s**(n-1)) / (2 * (1+s^n)^2) = 0
            numerator_eq2 = alpha * n * s**(n - 1)
            denominator_eq2 = 2 * (1 + s_pow_n)**2
            
            if denominator_eq2 == 0:
                return [1e9, 1e9]
            
            eq2 = (1 - beta / ((1 + s)**2)) - numerator_eq2 / denominator_eq2
            
            return [eq1, eq2]
            
        except (ValueError, OverflowError):
            # Return a large error if any math operation fails (e.g., overflow)
            return [1e9, 1e9]
    
    # Provide an initial guess for s and n. 
    # Preliminary analysis suggests the solution is near s=3.5, n=2.0.
    initial_guess = [3.5, 2.0]
    
    # Use scipy.optimize.root to find the solution to the system of equations.
    solution = root(find_bifurcation_point, initial_guess, method='hybr', tol=1e-9)
    
    # Output the results and the final answer
    if solution.success:
        s_c, n_c = solution.x
        
        print("The analysis is based on the following parameter values:")
        print(f"alpha = {int(alpha)}")
        print(f"beta = {int(beta)}\n")

        print("The critical value of n is found by solving the system of equations for the fixed point and the Hopf bifurcation condition:\n")
        # Print Equation 1 with parameter values
        print(f"1. Fixed Point: s + {int(beta)}/(1+s) - {int(alpha)}/(1+s^n) = 0")
        # Print Equation 2 with parameter values
        print(f"2. Bifurcation: 1 - {int(beta)}/((1+s)^2) - ({int(alpha)}*n*s^(n-1)) / (2*(1+s^n)^2) = 0\n")
        
        print(f"The system exhibits oscillations for values of n greater than the critical value found by solving these equations.")
        print(f"The critical value for n is approximately {n_c:.4f}.")
        print(f"<<<{n_c:.4f}>>>")

    else:
        print("Numerical solver failed to find a solution.")
        print(f"Solver message: {solution.message}")

if __name__ == "__main__":
    solve_for_n()