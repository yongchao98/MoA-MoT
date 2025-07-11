import numpy as np
from scipy.optimize import curve_fit

def solve():
    """
    This function finds the maximally parsimonious model for the given data,
    reports the final equation, and prints the estimated parameter values.
    """
    # Step 1: Define the dataset provided by the user
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # Step 2: Define the model function. Based on analysis, a square root
    # function y = a * sqrt(x) + b provides the most parsimonious and optimal fit.
    def sqrt_model(x, a, b):
        """A square root model of the form y = a*sqrt(x) + b."""
        return a * np.sqrt(x) + b

    # Step 3: Fit the model to the data to find the optimal parameters
    # 'params' will contain the optimal values for 'a' and 'b'.
    try:
        params, _ = curve_fit(sqrt_model, x_data, y_data)
        a_opt, b_opt = params
    except RuntimeError:
        print("Error: Optimal parameters not found. Could not fit the model.")
        return

    # Step 4: Format the parameters to 3 significant digits.
    # The .3g format specifier is ideal for this.
    a_str = f"{a_opt:.3g}"
    b_str = f"{b_opt:.3g}"
    
    # Step 5: Print the results in the required format.
    print("The maximally parsimonious model that optimally predicts y is of the form: y = a * sqrt(x) + b")
    
    # Handle the sign of the 'b' parameter for clean printing of the equation.
    if b_opt < 0:
        # Use abs() to remove the negative sign from b_str, as it's already included.
        final_equation = f"y = {a_str} * sqrt(x) - {abs(float(b_str))}"
    else:
        final_equation = f"y = {a_str} * sqrt(x) + {b_str}"

    print(f"\nFinal Equation: {final_equation}")

    print("\nParameter Values:")
    print(f"a = {a_str}")
    print(f"b = {b_str}")

solve()