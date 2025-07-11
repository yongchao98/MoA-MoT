import numpy as np
from scipy.optimize import minimize

def run_demonstration():
    """
    This function demonstrates that the statement 'Any strictly convex function has a
    unique global minimizer' is false using a counterexample.
    """

    # We consider the function f(x) = e^x.
    # Its second derivative is f''(x) = e^x, which is > 0 for all x.
    # A function with a strictly positive second derivative is strictly convex.
    def strictly_convex_function(x):
        return np.exp(x[0])

    print("Analyzing the statement: 'Any strictly convex function has a unique global minimizer'.")
    print("This statement is FALSE.\n")
    print("Let's use f(x) = e^x as a counterexample.")
    print("This function is strictly convex because its second derivative, e^x, is always positive.")
    print("However, this function has no global minimum. As x -> -infinity, f(x) -> 0, but it never reaches a minimum value.\n")
    print("We will now use a numerical optimizer to try and find the minimum.")
    
    # An initial starting point for the optimization
    initial_guess = np.array([10.0])
    
    # Attempt to minimize the function
    result = minimize(
        fun=strictly_convex_function,
        x0=initial_guess,
        method='BFGS' # A common optimization algorithm
    )

    print("\n--- Scipy Optimizer Result ---")
    print(f"Optimization Successful: {result.success}")
    print(f"Message: {result.message}")
    print("Note: A 'Desired error not necessarily achieved' message indicates failure to converge to a minimum.")
    print(f"Final x-value found: {result.x[0]}")
    print(f"Function value at final x: {result.fun}")
    print("----------------------------\n")

    print("The optimizer moved x to a large negative number and the function value approaches zero.")
    print("This demonstrates that the function does not have a finite global minimum.")
    print("Therefore, the statement that *any* strictly convex function has a unique global minimizer is false.")

run_demonstration()