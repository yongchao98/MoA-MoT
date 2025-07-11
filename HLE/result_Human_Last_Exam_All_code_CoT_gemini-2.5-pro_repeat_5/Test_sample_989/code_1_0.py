import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def demonstrate_convex_function_without_minimum():
    """
    This function demonstrates that a strictly convex function does not necessarily
    have a global minimizer, which addresses statement E.
    """
    # Define a strictly convex function: f(x) = e^x
    # The first derivative is e^x and the second derivative is e^x.
    # Since the second derivative is always > 0, the function is strictly convex.
    f = lambda x: np.exp(x)

    print("Analyzing statement E: 'Any strictly convex function has a unique global minimizer'.")
    print("Let's test this with the strictly convex function f(x) = e^x.")
    print("This function's value approaches 0 as x -> -inf, but it never reaches a minimum value.")
    print("\nAttempting to find the minimum using scipy.optimize.minimize...")

    # We start the search for the minimum at an arbitrary point, e.g., x=0.
    initial_guess = 0
    
    # The 'BFGS' method is a standard quasi-Newton method for unconstrained optimization.
    result = minimize(f, initial_guess, method='BFGS')

    print("\n--- Optimization Result ---")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Final function value (fun): {result.fun}")
    print(f"Final x value (x): {result.x[0]}")
    print("---------------------------\n")

    if not result.success:
        print("Conclusion: The optimization failed to converge to a minimum.")
        print("This supports the fact that f(x) = e^x, despite being strictly convex, has no global minimum.")
        print("Therefore, statement E is false.")
    else:
        # This case is unlikely but included for completeness.
        print("Conclusion: The optimizer claims to have found a minimum, but this is likely a numerical artifact.")
        print("Theoretically, f(x)=e^x has no minimum. Therefore, statement E is false.")

    # --- Plotting for visual confirmation ---
    x_vals = np.linspace(-10, 5, 400)
    y_vals = f(x_vals)

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, label='f(x) = e^x')
    plt.title('Plot of a Strictly Convex Function without a Global Minimum')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(True)
    plt.legend()
    plt.axhline(0, color='red', linestyle='--', label='Infimum (y=0)')
    plt.text(-9.5, 20, "This function is strictly convex (f''(x) = e^x > 0),\nbut it never reaches a minimum value.", fontsize=12)
    plt.text(-9.5, 5, "It only approaches an infimum of 0 as x approaches -infinity.", fontsize=12)
    
    # The plot is for visual understanding and will be displayed if run in an environment that supports it.
    # To make this script runnable in any environment, we won't call plt.show()
    # but the creation of the plot is part of the demonstration.
    print("\nA plot has been generated to visualize that f(x) = e^x has no minimum.")

demonstrate_convex_function_without_minimum()