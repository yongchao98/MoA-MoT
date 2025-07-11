import numpy as np

def solve_model():
    """
    Constructs a parsimonious model for the given dataset and prints the final equation.
    """
    # 1. Represent the Data
    # The data consists of 8 points from a 2^3 factorial design.
    x_data = np.array([
        [-1, -1, -1],
        [ 1, -1, -1],
        [-1,  1, -1],
        [ 1,  1, -1],
        [-1, -1,  1],
        [ 1, -1,  1],
        [-1,  1,  1],
        [ 1,  1,  1]
    ])
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])

    # Extract individual factors for clarity
    x1 = x_data[:, 0]
    x2 = x_data[:, 1]
    x3 = x_data[:, 2]

    # 2. Full Model Estimation
    # Construct the design matrix X for the full model including all interactions.
    # Columns correspond to: Intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3
    X = np.c_[
        np.ones(y.shape[0]),
        x1,
        x2,
        x3,
        x1 * x2,
        x1 * x3,
        x2 * x3,
        x1 * x2 * x3
    ]
    
    # 3. Calculate Coefficients
    # Using the least squares formula for an orthogonal design.
    n = y.shape[0]
    beta = (1/n) * X.T @ y

    # 4. Model Selection
    # Identify significant effects by comparing coefficient magnitudes.
    # The calculated coefficients are approximately:
    # b0=51.4, b1=20.8, b2=-0.0, b3=-12.8, b12=0.8, b13=0.1, b23=0.3, b123=2.7
    # The coefficients for x1 (b1), x3 (b3), and the three-way interaction x1*x2*x3 (b123)
    # are significantly larger in magnitude than the others (b2, b12, b13, b23).
    # Therefore, the parsimonious model will include the intercept, x1, x3, and x1*x2*x3.
    
    model_coeffs = {
        "intercept": beta[0],
        "x_1": beta[1],
        "x_3": beta[3],
        "x_1 * x_2 * x_3": beta[7]
    }

    # 5. Final Equation
    # Build the equation string from the selected model coefficients.
    
    # Start with the intercept
    equation = f"y = {model_coeffs['intercept']:.1f}"

    # Add the term for x1
    b1_val = model_coeffs['x_1']
    equation += f" {'+' if b1_val >= 0 else '-'} {abs(b1_val):.1f} * x_1"
    
    # Add the term for x3
    b3_val = model_coeffs['x_3']
    equation += f" {'+' if b3_val >= 0 else '-'} {abs(b3_val):.1f} * x_3"
    
    # Add the term for x1*x2*x3
    b123_val = model_coeffs['x_1 * x_2 * x_3']
    equation += f" {'+' if b123_val >= 0 else '-'} {abs(b123_val):.1f} * x_1 * x_2 * x_3"

    print(equation)

solve_model()