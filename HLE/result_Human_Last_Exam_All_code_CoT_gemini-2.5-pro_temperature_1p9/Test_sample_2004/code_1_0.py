import numpy as np

def solve():
    """
    This function devises a maximally parsimonious model for optimally predicting y
    from x using the provided dataset.
    """
    # The 25 observations of x and y
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    def r_squared(y_true, y_pred):
        """Calculates the R-squared value for a model's predictions."""
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - (ss_res / ss_tot)

    # --- Model Comparison ---
    # Model 1: Linear (y = ax + b), 2 parameters
    params_lin = np.polyfit(x_data, y_data, 1)
    r2_lin = r_squared(y_data, np.polyval(params_lin, x_data))

    # Model 2: Quadratic (y = ax^2 + bx + c), 3 parameters
    params_quad = np.polyfit(x_data, y_data, 2)
    r2_quad = r_squared(y_data, np.polyval(params_quad, x_data))

    # Model 3: Square Root (y = a*sqrt(x) + b), 2 parameters
    x_sqrt = np.sqrt(x_data)
    params_sqrt = np.polyfit(x_sqrt, y_data, 1)
    r2_sqrt = r_squared(y_data, np.polyval(params_sqrt, x_sqrt))

    # --- Model Selection ---
    # The square root model has the highest R-squared (0.958) while having only 2 parameters,
    # making it the most parsimonious choice.
    # R-squared values for comparison:
    # Linear: 0.866, Quadratic: 0.902, Square Root: 0.958

    # --- Final Model ---
    a, b = params_sqrt
    a_formatted = f"{a:.3g}"
    
    print("The maximally parsimonious model was found to be a square root function of the form: y = a*sqrt(x) + b")
    print("\nThe final equation with parameter estimates rounded to 3 significant digits is:")
    
    # Construct and print the final equation string with correct signs
    if b < 0:
        b_formatted = f"{abs(b):.3g}"
        print(f"y = {a_formatted} * sqrt(x) - {b_formatted}")
    else:
        b_formatted = f"{b:.3g}"
        print(f"y = {a_formatted} * sqrt(x) + {b_formatted}")

solve()
<<<y = 0.887 * sqrt(x) - 0.565>>>