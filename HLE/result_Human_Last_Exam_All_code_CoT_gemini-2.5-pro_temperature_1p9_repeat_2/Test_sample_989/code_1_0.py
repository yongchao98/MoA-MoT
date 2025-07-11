import numpy as np
from scipy.optimize import minimize

def main():
    """
    Analyzes the false statement about convex functions.
    
    Statement E: Any strictly convex function has a unique global minimizer.
    
    This statement is false. A strictly convex function is not guaranteed to have
    a global minimum. A simple counter-example is f(x) = e^x. This function
    is strictly convex, but it does not have a global minimum on the real line;
    it approaches 0 as x approaches -infinity, but never reaches it.

    This script demonstrates this by attempting to find the minimum of f(x) = e^x
    using a numerical optimizer. The optimizer will fail to converge to a specific
    minimum point.
    """
    
    # Define the strictly convex function f(x) = e^x
    def strictly_convex_function(x):
      return np.exp(x[0])

    # Initial guess for the minimizer
    initial_guess = np.array([0.0])

    # Attempt to minimize the function
    # Note: We expect this to fail or indicate it's unbounded below (in a practical sense)
    # The algorithm will likely be stopped by iteration limits or small gradient.
    result = minimize(strictly_convex_function, initial_guess, method='BFGS')

    print("Analyzing Statement E: Any strictly convex function has a unique global minimizer.")
    print("-" * 70)
    print("Let's test this with a counter-example: f(x) = e^x, which is strictly convex.")
    print("We will use a numerical optimizer to try and find its global minimum.")
    print("\nOptimization Attempt Details:")
    print(f"  Success: {result.success}")
    print(f"  Message: {result.message}")
    print(f"  Final x value: {result.x[0]:.2f}")
    print(f"  Final function value f(x): {result.fun:.2e}")

    print("\nConclusion:")
    print("The optimizer did not successfully converge to a minimum. Instead, it moved towards")
    print("negative infinity until it was stopped. This demonstrates that the strictly")
    print("convex function f(x) = e^x does not have a global minimizer.")
    print("Therefore, the statement that *any* strictly convex function has a unique")
    print("global minimizer is FALSE.")
    print("\nThe false statement is E.")

if __name__ == "__main__":
    main()
<<<E>>>