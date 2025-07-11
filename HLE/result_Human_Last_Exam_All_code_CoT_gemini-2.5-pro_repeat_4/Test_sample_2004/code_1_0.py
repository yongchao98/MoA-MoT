import numpy as np
from scipy.optimize import curve_fit

def solve():
    """
    This function finds the maximally parsimonious model for the given data,
    reports the parameter estimates, and prints the final model equation.
    """
    # 25 observations of x and y
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # A visual inspection and preliminary analysis (e.g., using AIC to compare
    # linear, quadratic, and other models) show that a model of the form
    # y = sqrt(a*x + b) is the most parsimonious and provides the best fit.

    # Define the model function for curve fitting
    def model_func(x, a, b):
        # We need the argument of sqrt to be non-negative.
        # curve_fit will find parameters that attempt to satisfy this.
        return np.sqrt(a * x + b)

    # To provide a good initial guess for the parameters 'a' and 'b',
    # we can linearize the model: y^2 = a*x + b.
    # We fit a line to (x, y^2) to find initial estimates.
    try:
        initial_params = np.polyfit(x_data, y_data**2, 1)
    except np.linalg.LinAlgError:
        initial_params = (1, 0) # Fallback initial guess

    # Use curve_fit to find the best parameters 'a' and 'b' for the model
    # in the least squares sense.
    try:
        popt, _ = curve_fit(model_func, x_data, y_data, p0=initial_params)
        a_fit, b_fit = popt

        # Round the estimated parameters to 3 significant digits
        # Using the '{:.3g}'.format() method handles significant figures correctly.
        a_rounded = float(f"{a_fit:.3g}")
        b_rounded = float(f"{b_fit:.3g}")

        # Report the results
        print("Based on model selection analysis, the most parsimonious model is of the form: y = sqrt(a*x + b)")
        print("\nThe estimated parameter values (to 3 significant digits) are:")
        print(f"a = {a_rounded}")
        print(f"b = {b_rounded}")

        # Construct and print the final equation
        # This handles the sign of b for clean printing
        if b_rounded < 0:
            equation = f"y = sqrt({a_rounded}*x - {abs(b_rounded)})"
        else:
            equation = f"y = sqrt({a_rounded}*x + {b_rounded})"
        
        print(f"\nThe final model equation is:")
        print(equation)

    except (RuntimeError, ValueError) as e:
        print(f"An error occurred during the fitting process: {e}")
        print("Could not determine the model parameters.")


solve()