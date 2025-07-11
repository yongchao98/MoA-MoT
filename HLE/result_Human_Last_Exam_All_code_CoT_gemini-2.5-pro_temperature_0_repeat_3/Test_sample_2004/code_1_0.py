import numpy as np

def solve():
    """
    This function finds the maximally parsimonious model for the given data,
    calculates its parameters, and prints the final equation.
    """
    # The 25 observations of x and y
    x_data = np.array([
        5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6,
        5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45
    ])
    y_data = np.array([
        1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166,
        1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112,
        1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123
    ])

    # Based on analysis, the square root model y = a*sqrt(x) + b is the best fit.
    # To find the parameters 'a' and 'b', we perform a linear regression
    # on y vs. sqrt(x).
    x_transformed = np.sqrt(x_data)

    # np.polyfit with degree 1 performs a linear least-squares regression.
    # It returns the coefficients [a, b] for the equation y = a*x_transformed + b.
    params = np.polyfit(x_transformed, y_data, 1)
    a = params[0]
    b = params[1]

    # Print the final equation with parameters formatted to 3 significant digits.
    # We handle the sign of the intercept 'b' for clean output.
    print("The most parsimonious model that optimally predicts y is:")
    if b >= 0:
        # Format parameters to 3 significant digits using the '.3g' specifier
        print(f"y = {a:.3g} * sqrt(x) + {b:.3g}")
    else:
        # If b is negative, print a minus sign instead of a plus.
        print(f"y = {a:.3g} * sqrt(x) - {abs(b):.3g}")

solve()