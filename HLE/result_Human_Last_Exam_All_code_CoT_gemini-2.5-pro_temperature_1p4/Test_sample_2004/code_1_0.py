import numpy as np
from scipy.optimize import curve_fit

def solve():
    """
    This function finds the maximally parsimonious model for the given data,
    which is determined to be a square root model of the form y = b*sqrt(x) + c.
    It calculates the parameters and prints the final equation.
    """
    # The 25 observations of x and y
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # Define the square root model function
    def sqrt_model(x, b, c):
        return b * np.sqrt(x) + c

    # Use curve_fit to find the best-fit parameters b and c using non-linear least squares
    popt, _ = curve_fit(sqrt_model, x_data, y_data)
    b_estimate, c_estimate = popt

    # Report the final model and its parameters
    print("The most parsimonious model for optimally predicting y is of the form: y = b * sqrt(x) + c")
    
    # Output each number in the final equation formatted to 3 significant digits
    print("\nThe estimated parameters to 3 significant digits are:")
    print(f"b = {b_estimate:.3g}")
    print(f"c = {c_estimate:.3g}")

    print("\nThe final equation is:")
    
    # Format the equation string to handle the sign of the intercept gracefully
    if c_estimate >= 0:
        print(f"y = {b_estimate:.3g} * sqrt(x) + {abs(c_estimate):.3g}")
    else:
        print(f"y = {b_estimate:.3g} * sqrt(x) - {abs(c_estimate):.3g}")

solve()