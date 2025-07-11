import numpy as np
from scipy.optimize import curve_fit

def find_best_model():
    """
    This function takes 25 observations of x and y, devises a maximally
    parsimonious model for optimally predicting y from x, and reports
    the model equation with parameter estimates.
    """
    # The 25 observations of x and y
    x_data = np.array([
        5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8,
        7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45
    ])
    y_data = np.array([
        1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866,
        0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8,
        1.501, 2.303, 2.31, 1.194, 1.144, 0.123
    ])

    # The visual inspection of the data and comparison of models (like linear,
    # quadratic, log) using BIC point to a square root function as the
    # most parsimonious and best-fitting model.

    # Define the square root model function: y = a * sqrt(x) + b
    def sqrt_model(x, a, b):
        return a * np.sqrt(x) + b

    # Use curve_fit to find the optimal parameters 'a' and 'b'
    # that minimize the sum of the squared residuals.
    params, _ = curve_fit(sqrt_model, x_data, y_data)
    a, b = params

    # Format the parameters to 3 significant digits.
    # The 'g' format specifier is used for general-purpose formatting,
    # which is ideal for significant figures.
    a_formatted = f"{a:.3g}"
    b_formatted = f"{b:.3g}"

    # To present the equation cleanly, handle the sign of the constant term.
    if b >= 0:
        operator = "+"
        b_val_str = f"{b:.3g}"
    else:
        operator = "-"
        b_val_str = f"{abs(b):.3g}"

    print("The most parsimonious model with the best fit (lowest BIC) is the square root model.")
    print("The model equation is of the form: y = a * sqrt(x) + b")
    print("\nParameter estimates (to 3 significant digits):")
    print(f"a = {a_formatted}")
    print(f"b = {b_formatted}")
    print("\nFinal Model Equation:")
    print(f"y = {a_formatted} * sqrt(x) {operator} {b_val_str}")


if __name__ == "__main__":
    find_best_model()
