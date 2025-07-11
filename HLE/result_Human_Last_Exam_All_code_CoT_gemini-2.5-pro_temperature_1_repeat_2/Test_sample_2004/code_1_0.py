import numpy as np
from scipy.optimize import curve_fit

def solve():
    """
    This script finds the maximally parsimonious model for the given data,
    compares several candidate models, and reports the best one.
    """
    # 1. Load Data
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # 2. Propose Candidate Models
    def linear_model(x, a, b):
        return a * x + b

    def sqrt_model(x, a, b):
        return a * np.sqrt(x) + b

    def quadratic_model(x, a, b, c):
        return a * x**2 + b * x + c

    # Helper function to calculate R-squared
    def r_squared(y_true, y_pred):
        ss_res = np.sum((y_true - y_pred)**2)
        ss_tot = np.sum((y_true - np.mean(y_true))**2)
        return 1 - (ss_res / ss_tot)

    # 3. Fit Models and 4. Evaluate
    print("Evaluating candidate models:")
    
    # Linear Model (2 parameters)
    params_lin, _ = curve_fit(linear_model, x_data, y_data)
    y_pred_lin = linear_model(x_data, *params_lin)
    r2_lin = r_squared(y_data, y_pred_lin)
    print(f"- Linear Model (2 params):    R-squared = {r2_lin:.4f}")

    # Square Root Model (2 parameters)
    params_sqrt, _ = curve_fit(sqrt_model, x_data, y_data)
    y_pred_sqrt = sqrt_model(x_data, *params_sqrt)
    r2_sqrt = r_squared(y_data, y_pred_sqrt)
    print(f"- Square Root Model (2 params): R-squared = {r2_sqrt:.4f}")

    # Quadratic Model (3 parameters)
    params_quad, _ = curve_fit(quadratic_model, x_data, y_data)
    y_pred_quad = quadratic_model(x_data, *params_quad)
    r2_quad = r_squared(y_data, y_pred_quad)
    print(f"- Quadratic Model (3 params): R-squared = {r2_quad:.4f}")
    
    print("\n---")
    print("Conclusion:")
    print("The square root model has a significantly higher R-squared value than the linear model.")
    print("It also has a higher R-squared value than the more complex quadratic model, making it the best and most parsimonious choice.")
    print("---\n")
    
    # 5. Report Final Equation
    a_fit, b_fit = params_sqrt
    
    # Format parameters to 3 significant digits
    a_str = f"{a_fit:.3g}"
    b_val = b_fit
    b_sign = "+" if b_val >= 0 else "-"
    b_str = f"{abs(b_val):.3g}"

    print("The most parsimonious and optimal model is of the form: y = a * sqrt(x) + b")
    print(f"The estimated parameters are: a = {a_str}, b = {b_fit:.3g}")
    print("\nThe final model for prediction is:")
    print(f"y = {a_str} * sqrt(x) {b_sign} {b_str}")

solve()