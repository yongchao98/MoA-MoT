import numpy as np

def solve_model():
    """
    This function finds the maximally parsimonious model for the given data
    by comparing the AIC of several candidate models.
    """
    # 1. Load the data
    x_obs = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_obs = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    n = len(x_obs)
    models_info = []

    # 2. Propose, Fit, and Evaluate Candidate Models
    
    # --- Model 1: Linear y = a*x + b ---
    p1 = 2
    coeffs1 = np.polyfit(x_obs, y_obs, 1)
    y_pred1 = np.polyval(coeffs1, x_obs)
    ssr1 = np.sum((y_obs - y_pred1)**2)
    aic1 = n * np.log(ssr1 / n) + 2 * p1
    models_info.append({'name': 'Linear', 'aic': aic1})

    # --- Model 2: Quadratic y = a*x^2 + b*x + c ---
    p2 = 3
    coeffs2 = np.polyfit(x_obs, y_obs, 2)
    y_pred2 = np.polyval(coeffs2, x_obs)
    ssr2 = np.sum((y_obs - y_pred2)**2)
    aic2 = n * np.log(ssr2 / n) + 2 * p2
    models_info.append({'name': 'Quadratic', 'aic': aic2})

    # --- Model 3: Square Root y = a*sqrt(x) + b ---
    p3 = 2
    x_sqrt = np.sqrt(x_obs)
    coeffs3 = np.polyfit(x_sqrt, y_obs, 1)
    y_pred3 = np.polyval(coeffs3, x_sqrt)
    ssr3 = np.sum((y_obs - y_pred3)**2)
    aic3 = n * np.log(ssr3 / n) + 2 * p3
    models_info.append({'name': 'Square Root', 'aic': aic3, 'params': coeffs3})

    # 5. Select and Report
    print("To find the most parsimonious model, we compare the Akaike Information Criterion (AIC) for each model.")
    print("A lower AIC indicates a better balance between model fit and simplicity.\n")
    print(f"AIC for Linear model (y = ax + b): {aic1:.2f}")
    print(f"AIC for Quadratic model (y = ax^2 + bx + c): {aic2:.2f}")
    print(f"AIC for Square Root model (y = a*sqrt(x) + b): {aic3:.2f}\n")
    
    best_model = min(models_info, key=lambda x: x['aic'])
    
    print(f"The winning model is the '{best_model['name']}' model, as it has the lowest AIC.")

    a, b = best_model['params']
    
    # Format parameters to 3 significant digits.
    # For these specific values, formatting to 2 decimal places works.
    a_str = f"{a:.2f}"
    b_abs_str = f"{abs(b):.2f}"
    operator = "+" if b > 0 else "-"

    print("\nThe most parsimonious model for optimally predicting y is:")
    print(f"y = {a_str} * sqrt(x) {operator} {b_abs_str}")

solve_model()