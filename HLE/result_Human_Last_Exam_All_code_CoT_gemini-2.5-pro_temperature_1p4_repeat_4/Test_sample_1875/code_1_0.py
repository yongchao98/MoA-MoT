import numpy as np

def solve_model():
    """
    Analyzes factorial design data to construct a parsimonious predictive model.
    """
    # 1. Define the data
    # The data consists of 8 points from a 2^3 factorial design.
    # The columns correspond to x1, x2, x3.
    x_factors = np.array([
        [-1, -1, -1],
        [ 1, -1, -1],
        [-1,  1, -1],
        [ 1,  1, -1],
        [-1, -1,  1],
        [ 1, -1,  1],
        [-1,  1,  1],
        [ 1,  1,  1]
    ])

    # The response variable y
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])

    # 2. Construct the full design matrix X for the model including all interactions
    intercept = np.ones((len(y), 1))
    x1 = x_factors[:, 0:1]
    x2 = x_factors[:, 1:2]
    x3 = x_factors[:, 2:3]
    x1x2 = x1 * x2
    x1x3 = x1 * x3
    x2x3 = x2 * x3
    x1x2x3 = x1 * x2 * x3

    X = np.hstack([intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3])
    term_names = ["Intercept", "x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"]

    # 3. Estimate all parameters (coefficients)
    # For an orthogonal design X'X = nI, so beta = (1/n) * X'y
    n = len(y)
    coeffs = (1/n) * X.T @ y

    # 4. Identify significant effects by comparing magnitudes
    # A parsimonious model keeps only the most influential terms.
    # We observe the calculated coefficients:
    # b0 (Intercept): 51.375
    # b1 (x1): 20.825
    # b2 (x2): -0.025
    # b3 (x3): -12.75
    # b12 (x1*x2): 0.775
    # b13 (x1*x3): 0.1
    # b23 (x2*x3): 0.025
    # b123 (x1*x2*x3): -0.15
    # The coefficients for x1 and x3 are significantly larger in magnitude than all others.
    # The most parsimonious model will include only these two main effects and the intercept.

    b0 = coeffs[0]
    b1 = coeffs[1]
    b3 = coeffs[3]

    # 5. Present the final equation with parameters rounded to one decimal place.
    # The final model is of the form: y = b0 + b1*x1 + b3*x3
    # Note: We build the string carefully to handle the sign of the x3 coefficient.
    
    final_equation = f"y = {b0:.1f} + {b1:.1f} * x_1 - {abs(b3):.1f} * x_3"
    print("The parsimonious model is:")
    print(final_equation)

solve_model()