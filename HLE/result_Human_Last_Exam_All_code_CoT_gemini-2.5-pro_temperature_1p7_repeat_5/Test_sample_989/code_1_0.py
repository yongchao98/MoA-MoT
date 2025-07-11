import numpy as np
from scipy.optimize import minimize
import warnings

# Suppress the optimization warning for a cleaner output
warnings.filterwarnings("ignore", message="Desired error not necessarily achieved due to precision loss.")

def run_demonstration():
    """
    This function demonstrates that the statement "Any strictly convex function
    has a unique global minimizer" is false by providing a counterexample.
    """

    # Define the strictly convex function f(x) = e^x
    # It is strictly convex because its second derivative, e^x, is always positive.
    def strictly_convex_function(x):
        return np.exp(x[0])

    print("Analyzing the statement: 'Any strictly convex function has a unique global minimizer.'")
    print("Let's test this with the function f(x) = e^x, which is strictly convex.")
    print("The infimum (greatest lower bound) of this function is 0, but this value is never achieved for any finite x.")
    print("Therefore, the function has no global minimizer.\n")
    print("We will now use a numerical optimizer from SciPy to try and find the minimum.")
    print("We expect the optimizer to fail or to search endlessly towards negative infinity.\n")

    # An initial starting point for the optimization
    start_point = np.array([10.0])

    # Attempt to find the minimum of the function.
    # The 'BFGS' algorithm is a standard method for unconstrained minimization.
    result = minimize(
        fun=strictly_convex_function,
        x0=start_point,
        method='BFGS'
    )

    print("--- Optimizer Result ---")
    print(f"Success: {result.success}")
    print(f"Status Message: {result.message}")
    print(f"Final x value: {result.x[0]}")
    print(f"Final function value (e^x): {result.fun}")
    print("------------------------\n")
    print("Conclusion:")
    print("The optimizer did not successfully converge. It stopped after heading towards negative infinity,")
    print("achieving a very small function value, but this is not a true minimum.")
    print("This demonstrates that a strictly convex function does not necessarily have a minimizer.")
    print("Therefore, the statement is FALSE.")


if __name__ == '__main__':
    run_demonstration()
