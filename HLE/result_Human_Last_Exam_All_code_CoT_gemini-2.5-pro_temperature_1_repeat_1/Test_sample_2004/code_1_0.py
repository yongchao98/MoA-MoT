import numpy as np
from scipy.optimize import curve_fit

def solve():
    """
    This function finds the maximally parsimonious model for the given data,
    estimates the parameters, and prints the final equation.
    """
    # The 25 observations of x and y
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # Define the square root model function, which was determined to be the best fit.
    # The function takes the independent variable 'x' and the parameters 'a' and 'b' to be fitted.
    def sqrt_model(x, a, b):
        return a * np.sqrt(x) + b

    # Use scipy's curve_fit to find the optimal parameters 'a' and 'b'
    # that minimize the sum of the squared residuals.
    # popt contains the optimal parameters [a, b].
    popt, _ = curve_fit(sqrt_model, x_data, y_data)
    a_fit, b_fit = popt

    # Format the estimated parameters to 3 significant digits.
    # The format specifier '{:.3g}' is used for this purpose.
    a_formatted = '{:.3g}'.format(a_fit)
    b_formatted = '{:.3g}'.format(b_fit)
    
    # Print the results, including the final model equation.
    print("The best-fit model is a square root function of the form: y = a * sqrt(x) + b")
    print("\nEstimated parameter values (to 3 significant digits):")
    print(f"a = {a_formatted}")
    print(f"b = {b_formatted}")

    print("\nThe final prediction equation is:")
    # We check if 'b' is negative to correctly format the sign in the equation.
    if b_fit < 0:
        # Use abs() to remove the negative sign from 'b' as it's already included.
        print(f"y = {a_formatted} * sqrt(x) - {abs(float(b_formatted))}")
    else:
        print(f"y = {a_formatted} * sqrt(x) + {b_formatted}")

# Execute the function to solve the problem
solve()