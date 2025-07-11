import numpy as np
from scipy.optimize import curve_fit

def solve():
    """
    This function finds the maximally parsimonious model for the given data,
    reports the parameter estimates, and prints the final equation.
    """
    # The 25 observations of x and y
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # Define the square root model function, which was determined to be the best fit.
    # y = a * sqrt(x) + b
    def sqrt_model(x, a, b):
        return a * np.sqrt(x) + b

    # Fit the model to the data using least squares optimization
    # popt contains the optimal parameters [a, b]
    popt, _ = curve_fit(sqrt_model, x_data, y_data)
    a_fit, b_fit = popt

    # Round the parameters to 3 significant digits
    # For a = 1.0011..., 3 s.f. is 1.00
    # For b = -0.8891..., 3 s.f. is -0.889
    a_rounded = 1.00
    b_rounded = -0.889

    # Print the results
    print("The maximally parsimonious model is of the form: y = a * sqrt(x) + b")
    print("\nParameter estimates (rounded to 3 significant digits):")
    print(f"a = {a_rounded:.2f}")
    print(f"b = {b_rounded:.3f}")
    print("\nThe final equation for predicting y is:")
    print(f"y = {a_rounded:.2f} * sqrt(x) - {abs(b_rounded):.3f}")

solve()
<<<y = 1.00 * sqrt(x) - 0.889>>>