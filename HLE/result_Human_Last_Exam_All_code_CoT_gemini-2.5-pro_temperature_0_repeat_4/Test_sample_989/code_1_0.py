import numpy as np
from scipy.optimize import minimize

# The user wants to identify the false statement about SVMs.
# After analysis, statement E is identified as false.
# Statement E: "Any strictly convex function has a unique global minimizer"
# This script will provide a counterexample to prove that statement E is false.

def demonstrate_falsehood_of_statement_e():
    """
    This function demonstrates that statement E is false by using a counterexample.
    The function f(x) = e^x is strictly convex, but it does not have a global minimum.
    """
    print("--- Testing Statement E: 'Any strictly convex function has a unique global minimizer' ---")
    print("\nWe will use the function f(x) = e^x as a counterexample.")
    print("A function is strictly convex if its second derivative is always positive.")
    print("For f(x) = e^x, the second derivative is f''(x) = e^x, which is > 0 for all real numbers x.")
    print("Therefore, f(x) = e^x is a strictly convex function.")

    # Define the strictly convex function
    def convex_function(x):
        return np.exp(x[0])

    print("\nNow, we will use a numerical optimizer to try and find a minimum for this function.")
    print("If statement E were true, the optimizer should find a unique global minimum.")

    # An initial guess for the optimizer
    initial_guess = [0]

    # Run the optimizer. We expect it to fail to find a minimum.
    # The default 'BFGS' algorithm might struggle and warn.
    result = minimize(convex_function, initial_guess, method='BFGS')

    print(f"\n--- Optimizer Result ---")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Final function value: {result.fun}")
    print(f"Final x value: {result.x}")

    print("\n--- Analysis ---")
    if not result.success:
        print("As expected, the optimizer failed to converge to a minimum.")
        print("The function f(x) = e^x continuously decreases as x approaches negative infinity.")
        print("It approaches a lower bound (infimum) of 0, but never reaches it.")
        print("There is no specific x value that produces a minimum value for the function.")
    else:
        # This case is unlikely but depends on the optimizer's stopping criteria.
        print("The optimizer reported success, but this is likely due to reaching a tolerance limit.")
        print("The function value is very small, but an even smaller value can always be found by choosing a more negative x.")

    print("\n--- Conclusion ---")
    print("The existence of a strictly convex function (e^x) without a global minimizer proves that")
    print("Statement E is FALSE.")

if __name__ == '__main__':
    demonstrate_falsehood_of_statement_e()