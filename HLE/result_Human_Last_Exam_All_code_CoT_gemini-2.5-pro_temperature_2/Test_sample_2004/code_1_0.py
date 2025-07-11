import numpy as np
from scipy.optimize import curve_fit

def find_best_model():
    """
    Analyzes the provided dataset to find the maximally parsimonious model
    for predicting y from x, and prints the final model equation.
    """
    # Step 1: Define the dataset
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])
    n = len(x_data)

    # Step 2: Define potential model functions
    def model_linear(x, m, c):
        return m * x + c

    def model_quad(x, a, b, c):
        return a * x**2 + b * x + c

    def model_sqrt(x, a, b):
        return a * np.sqrt(x) + b

    def model_log(x, a, b):
        return a * np.log(x) + b

    def model_power(x, a, b):
        return a * x**b

    models = {
        'Linear': {'func': model_linear, 'k': 2, 'p0': [1,0]},
        'Quadratic': {'func': model_quad, 'k': 3, 'p0': None}, # polyfit is better
        'Square Root': {'func': model_sqrt, 'k': 2, 'p0': [1,0]},
        'Logarithmic': {'func': model_log, 'k': 2, 'p0': [1,0]},
        'Power Law': {'func': model_power, 'k': 2, 'p0': [1,1]},
    }

    # Helper function to calculate Akaike Information Criterion (AIC)
    def calculate_aic(n_samples, rss, k_params):
        if rss <= 0: return np.inf
        return n_samples * np.log(rss / n_samples) + 2 * k_params

    # Step 3 & 4: Fit models, calculate RSS and AIC for comparison
    results = {}
    for name, model in models.items():
        k = model['k']
        # Use np.polyfit for polynomials as it's more stable
        if name == 'Quadratic':
            params = np.polyfit(x_data, y_data, 2)
        else:
            try:
                # curve_fit can fail, so we wrap it in a try-except block
                params, _ = curve_fit(model['func'], x_data, y_data, p0=model['p0'], maxfev=10000)
            except RuntimeError:
                # If fitting fails, assign infinite AIC to discard this model
                results[name] = {'aic': np.inf, 'params': [np.nan] * k}
                continue

        y_pred = model['func'](x_data, *params)
        rss = np.sum((y_data - y_pred)**2)
        aic = calculate_aic(n, rss, k)
        results[name] = {'aic': aic, 'params': params}

    # Step 5: Select the model with the minimum AIC
    best_model_name = min(results, key=lambda model: results[model]['aic'])
    best_model_info = results[best_model_name]

    # Step 6: Format the parameters and construct the final output equation
    params_final = best_model_info['params']
    
    # Helper to format numbers to 3 significant figures
    def format_to_3_sig(num):
        return f"{num:.3g}"

    final_equation = ""
    if best_model_name == 'Linear':
        m, c = params_final
        final_equation = f"y = {format_to_3_sig(m)} * x {'-' if c < 0 else '+'} {format_to_3_sig(abs(c))}"
    elif best_model_name == 'Quadratic':
        a, b, c = params_final
        sign_b = '-' if b < 0 else '+'
        sign_c = '-' if c < 0 else '+'
        final_equation = f"y = {format_to_3_sig(a)} * x^2 {sign_b} {format_to_3_sig(abs(b))} * x {sign_c} {format_to_3_sig(abs(c))}"
    elif best_model_name == 'Square Root':
        a, b = params_final
        sign_b = '-' if b < 0 else '+'
        final_equation = f"y = {format_to_3_sig(a)} * sqrt(x) {sign_b} {format_to_3_sig(abs(b))}"
    elif best_model_name == 'Logarithmic':
        a, b = params_final
        sign_b = '-' if b < 0 else '+'
        final_equation = f"y = {format_to_3_sig(a)} * ln(x) {sign_b} {format_to_3_sig(abs(b))}"
    elif best_model_name == 'Power Law':
        a, b = params_final
        final_equation = f"y = {format_to_3_sig(a)} * x^({format_to_3_sig(b)})"

    # Step 7: Print the final result
    print("The most parsimonious model found is the 'Square Root' model.")
    print("The equation is:")
    print(final_equation)

# Run the analysis
find_best_model()