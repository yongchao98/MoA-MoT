import numpy as np

def solve_parsimonious_model():
    """
    Constructs and prints a parsimonious model to predict y from x1, x2, and x3.
    """
    # Define the provided data points
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

    x1 = data[:, 0]
    x2 = data[:, 1]
    x3 = data[:, 2]
    y = data[:, 3]

    # Construct the design matrix for the full model including all interactions.
    # The columns are: 1 (intercept), x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3
    X = np.c_[
        np.ones(len(data)),
        x1,
        x2,
        x3,
        x1 * x2,
        x1 * x3,
        x2 * x3,
        x1 * x2 * x3
    ]
    
    term_names = ["", "x_1", "x_2", "x_3", "x_1x_2", "x_1x_3", "x_2x_3", "x_1x_2x_3"]

    # Solve the linear system y = X*beta for beta using least squares.
    # beta = (X^T * X)^-1 * X^T * y
    coeffs = np.linalg.solve(X.T @ X, X.T @ y)

    # Build the parsimonious model equation string.
    # We start with the intercept.
    equation_str = f"y = {coeffs[0]:.1f}"

    # Add other terms if their coefficient's magnitude is significant.
    # A heuristic cutoff of 1.0 is used, as there is a clear gap in coefficient
    # magnitudes above and below this value.
    for i in range(1, len(coeffs)):
        if abs(coeffs[i]) > 1.0:
            coeff_val = coeffs[i]
            term_name = term_names[i]
            
            sign = "-" if coeff_val < 0 else "+"
            abs_coeff = abs(coeff_val)
            
            equation_str += f" {sign} {abs_coeff:.1f}{term_name}"
            
    # Print the final equation
    print(equation_str)

solve_parsimonious_model()
<<<y = 51.4 + 20.8x_1 - 12.6x_3 - 10.3x_2x_3 + 10.2x_1x_2x_3>>>