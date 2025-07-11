import numpy as np

def solve():
    """
    This function finds the maximally parsimonious model for the given x and y data,
    and prints the resulting equation with parameters rounded to 3 significant digits.
    """
    # The 25 observations of x and y
    x_data = np.array([
        5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6,
        7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45
    ])
    y_data = np.array([
        1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93,
        1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123
    ])

    # The chosen model is y = a * sqrt(x) + b, as it provides the best parsimonious fit.
    # We transform x to sqrt(x) to linearize the model.
    x_transformed = np.sqrt(x_data)

    # Perform a linear regression (polynomial fit of degree 1) on the transformed data.
    # This finds the optimal 'a' (slope) and 'b' (intercept) in the least-squares sense.
    params = np.polyfit(x_transformed, y_data, 1)
    a, b = params

    # Format the parameters to 3 significant digits.
    # Use an f-string with the 'g' format specifier for significant figures.
    a_formatted = f"{a:.3g}"
    
    # Format the equation string to handle the sign of the intercept gracefully.
    sign = '-' if b < 0 else '+'
    b_abs_formatted = f"{abs(b):.3g}"
    
    # Print the final equation, which represents the model.
    # The print statement outputs each number as requested.
    print(f"y = {a_formatted} * sqrt(x) {sign} {b_abs_formatted}")

solve()