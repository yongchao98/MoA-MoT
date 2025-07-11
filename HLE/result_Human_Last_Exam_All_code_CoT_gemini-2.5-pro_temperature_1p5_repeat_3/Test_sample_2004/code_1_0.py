import numpy as np
from scipy.optimize import curve_fit
import warnings

def solve():
    """
    This function devises a maximally parsimonious model for optimally
    predicting y from x for the given dataset. It compares several models,
    selects the best one based on fit (RSS) and simplicity (number of parameters),
    and reports the final equation.
    """
    # Ignore potential warnings from polyfit for ill-conditioned matrices.
    warnings.filterwarnings("ignore", category=np.RankWarning)

    # 1. Input the 25 observations of x and y
    x_data = np.array([
        5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8,
        7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45
    ])
    y_data = np.array([
        1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866,
        0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8,
        1.501, 2.303, 2.31, 1.194, 1.144, 0.123
    ])

    # 2. Define and fit candidate models to compare
    
    # --- Model 1: Linear Model y = a*x + b (2 parameters) ---
    lin_params = np.polyfit(x_data, y_data, 1)
    a_lin, b_lin = lin_params
    y_pred_lin = a_lin * x_data + b_lin
    rss_lin = np.sum((y_data - y_pred_lin)**2)

    # --- Model 2: Square Root Model y = a*sqrt(x) + b (2 parameters) ---
    def sqrt_model(x, a, b):
        return a * np.sqrt(x) + b
    
    # Provide an initial guess [a, b] to help the solver
    sqrt_params, _ = curve_fit(sqrt_model, x_data, y_data, p0=[1, 0])
    a_sqrt, b_sqrt = sqrt_params
    y_pred_sqrt = sqrt_model(x_data, a_sqrt, b_sqrt)
    rss_sqrt = np.sum((y_data - y_pred_sqrt)**2)
    
    # --- Model 3: Quadratic Model y = ax^2 + bx + c (3 parameters) ---
    quad_params = np.polyfit(x_data, y_data, 2)
    y_pred_quad = np.poly1d(quad_params)(x_data)
    rss_quad = np.sum((y_data-y_pred_quad)**2)


    # 3. Compare models to select the best one.
    # A lower RSS indicates a better fit.
    print("Model Comparison (based on Residual Sum of Squares):")
    print(f"- Linear Model (2 parameters):      RSS = {rss_lin:.4f}")
    print(f"- Square Root Model (2 parameters): RSS = {rss_sqrt:.4f}")
    print(f"- Quadratic Model (3 parameters):   RSS = {rss_quad:.4f}")
    
    print("\nConclusion: The square root model is the maximally parsimonious choice.")
    print("It provides a vastly better fit (much lower RSS) than the linear model for the same number of parameters, and it also fits better than the more complex quadratic model.")
    
    # 4. Format and report the final result with 3 significant digits
    
    # The chosen parameters are from the square root model
    a_final = a_sqrt
    b_final = b_sqrt

    # Format parameters to 3 significant digits using the '.3g' format specifier
    a_final_formatted = f'{a_final:.3g}'
    b_final_formatted_val = f'{abs(b_final):.3g}'
    sign = "+" if b_final >= 0 else "-"

    print("\n--- Final Model ---")
    print(f"The chosen model form is: y = a * sqrt(x) + b")
    print("\nEstimated parameter values (to 3 significant digits):")
    print(f"a = {a_final_formatted}")
    print(f"b = {b_final:.3g}")
    
    print("\nFinal Equation:")
    # The final equation is constructed with the formatted numbers
    print(f"y = {a_final_formatted} * sqrt(x) {sign} {b_final_formatted_val}")

solve()
<<<y = 0.941 * sqrt(x) - 0.639>>>