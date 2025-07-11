import numpy as np
from scipy.optimize import root_scalar

def verify_solution(z, k):
    """
    Verifies if z is a solution for a given branch k.
    The equation is z * i = i^z.
    i^z is calculated as exp(z * log(i)), where log(i) corresponds to branch k.
    """
    lhs = z * 1j
    
    # Calculate ln(i) for branch k
    log_i_k = 1j * (np.pi/2 + 2 * k * np.pi)
    
    # Calculate i^z for branch k
    rhs = np.exp(z * log_i_k)
    
    # Check if the sides are close enough
    is_solution = np.isclose(lhs, rhs)
    
    print(f"\nVerifying solution z = {z} for k = {k}:")
    print(f"  Equation: ({z}) * i = i^({z})")
    print(f"  Left side (z*i): {lhs}")
    print(f"  Right side (i^z): {rhs}")
    print(f"  Is it a solution for branch k={k}? {is_solution}")

def solve_and_print():
    """
    Finds and prints solutions to z*i = i^z.
    """
    print("Solving the complex equation: z * i = i^z")
    
    # --- Real Solutions ---
    print("\n--- Real Solutions ---")
    print("We found two real solutions: z = 1 and z = -1.")
    
    # z = 1
    # For z=1, we need 1*i = i^1. This holds for the principal branch k=0.
    verify_solution(1, k=0)

    # z = -1
    # For z=-1, we need -1*i = i^-1. This holds for the principal branch k=0.
    verify_solution(-1, k=0)
    
    # --- Pure Imaginary Solutions ---
    print("\n--- Pure Imaginary Solutions ---")
    print("An infinite series of pure imaginary solutions z = i*y exists for integers k <= -1.")
    print("y is the solution to the equation: -y = exp(-y * (pi/2 + 2*k*pi)).")
    print("Finding the first 5 solutions numerically:")
    
    for k in range(-1, -6, -1):
        theta_k = np.pi/2 + 2 * k * np.pi
        
        # Define the function whose root we want to find
        # f(y) = y + exp(-y * theta_k)
        func = lambda y: y + np.exp(-y * theta_k)
        
        # We know the solution is in the interval (-1, 0)
        try:
            sol = root_scalar(func, bracket=(-1, 0), method='brentq')
            if sol.converged:
                y_k = sol.root
                z_k = y_k * 1j
                verify_solution(z_k, k=k)
            else:
                print(f"\nSolver did not converge for k = {k}")
        except ValueError as e:
            print(f"\nCould not find a root for k = {k} in the given interval. {e}")

    print("\nNote: Infinitely more pure imaginary solutions exist for k < -5.")
    print("Other complex solutions (with non-zero real and imaginary parts) also exist but are harder to find.")

if __name__ == '__main__':
    solve_and_print()