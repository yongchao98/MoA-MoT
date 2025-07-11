import numpy as np

def solve_parsimonious_model():
    """
    This function constructs a parsimonious model for the given dataset by
    fitting a full model and then selecting the most significant terms.
    """
    # Step 1: Define the data points from the problem description.
    # The data is a 2-level, 3-factor full factorial design.
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])

    # Step 2: Build the design matrix X for the full model.
    # The model includes an intercept, main effects, 2-way interactions, and a 3-way interaction.
    intercept = np.ones(y.shape[0])
    x1x2 = x1 * x2
    x1x3 = x1 * x3
    x2x3 = x2 * x3
    x1x2x3 = x1 * x2 * x3

    # The order of columns is: intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3
    X = np.column_stack([intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3])

    # Step 3: Solve for the coefficients of the full model using least squares.
    # With 8 data points and 8 parameters, this provides an exact fit.
    coefficients, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

    # Step 4: To build a parsimonious model, we select the terms with the largest effect sizes (coefficients).
    # Based on the magnitudes, the three most significant effects are x1, x3, and the x2*x3 interaction.
    # We will build our final model with the intercept and these three terms.
    
    # Extract the coefficients for the selected terms and round to one decimal place.
    b0 = round(coefficients[0], 1)   # Intercept
    b1 = round(coefficients[1], 1)   # Coefficient for x1
    b3 = round(coefficients[3], 1)   # Coefficient for x3
    b23 = round(coefficients[6], 1) # Coefficient for x2*x3

    # Step 5: Construct the final model equation string.
    # The format will be y = f(x), showing each term and its estimated parameter.
    
    # Start with the intercept
    equation = f"y = {b0}"

    # Append the x1 term, handling the sign for proper formatting.
    if b1 >= 0:
        equation += f" + {b1}*x1"
    else:
        equation += f" - {abs(b1)}*x1"

    # Append the x3 term
    if b3 >= 0:
        equation += f" + {b3}*x3"
    else:
        equation += f" - {abs(b3)}*x3"

    # Append the x2*x3 interaction term
    if b23 >= 0:
        equation += f" + {b23}*x2*x3"
    else:
        equation += f" - {abs(b23)}*x2*x3"
        
    print(equation)

solve_parsimonious_model()