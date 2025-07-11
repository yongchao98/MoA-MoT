import numpy as np
from scipy.optimize import curve_fit

def devise_and_report_model():
    """
    This function identifies the maximally parsimonious model for the given data,
    calculates its parameters, and prints the final equation.
    """
    # 25 observations of x and y
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # Based on analysis, the square root model is the most parsimonious.
    # Model: y = a * sqrt(x) + b
    def sqrt_model(x, a, b):
        return a * np.sqrt(x) + b

    # Fit the model to the data using least squares optimization
    try:
        params, _ = curve_fit(sqrt_model, x_data, y_data)
        a, b = params
    except RuntimeError:
        print("Error: Could not find optimal parameters for the model.")
        return

    # Format parameters to 3 significant digits
    # The '.3g' format specifier handles this automatically.
    a_str = f"{a:.3g}"
    
    # Build the final equation string with clean signs
    if b >= 0:
        op = "+"
        b_str = f"{b:.3g}"
    else:
        op = "-"
        b_str = f"{abs(b):.3g}"

    print("The maximally parsimonious model was found to be a square root function of the form: y = a * sqrt(x) + b")
    print("\nEstimated parameter values (to 3 significant digits):")
    print(f"a = {a_str}")
    print(f"b = {b if b >= 0 else -float(b_str)}") # Print b with its sign

    print("\nThe final model for predicting y from x is:")
    print(f"y = {a_str} * sqrt(x) {op} {b_str}")

# Execute the function
devise_and_report_model()