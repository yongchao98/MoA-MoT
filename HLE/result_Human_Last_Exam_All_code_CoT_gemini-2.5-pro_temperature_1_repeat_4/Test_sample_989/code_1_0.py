import numpy as np
from scipy.optimize import minimize

def run_analysis():
    """
    Analyzes the statement: 'Any strictly convex function has a unique global minimizer.'
    """
    print("Analyzing Statement E: 'Any strictly convex function has a unique global minimizer.'")
    print("-" * 70)
    print("We will test this statement with a counterexample: the function f(x) = e^x.")
    print("A function is strictly convex if its second derivative is always positive.")
    print("For f(x) = e^x, the first derivative is f'(x) = e^x, and the second derivative is f''(x) = e^x.")
    print("Since f''(x) = e^x is > 0 for all real x, the function is strictly convex.")

    print("\nNow, let's try to find its global minimum using a numerical optimizer from SciPy.")
    
    # Define the strictly convex function
    def f(x):
        return np.exp(x)

    # Initial guess for the optimizer
    x0 = 0
    
    # Attempt to find the minimum
    # For a convex function, any local minimum found would be the global minimum.
    result = minimize(f, x0, method='BFGS')

    print("\n--- Optimizer Result ---")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Final x value: {result.x[0]:.2f}")
    print(f"Final function value (minimum found): {result.fun:.2e}")
    print("------------------------\n")

    print("The optimizer has stopped, but the result does not represent a true minimum.")
    print("The function f(x) = e^x approaches 0 as x approaches -infinity, but it never reaches 0.")
    print("It has an infimum (greatest lower bound) of 0, but no minimum value.")
    print("This demonstrates that a strictly convex function does not necessarily have a global minimizer.")
    print("\nConclusion: Statement E is false.")

if __name__ == '__main__':
    run_analysis()