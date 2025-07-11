import numpy as np

def create_parsimonious_model():
    """
    Analyzes the provided experimental data to construct a parsimonious
    predictive model for y.
    """
    # Data from the problem description
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

    y = data[:, 3]
    x1, x2, x3 = data[:, 0], data[:, 1], data[:, 2]

    # Construct the design matrix for the full model including all interactions
    X_full = np.c_[
        np.ones(len(x1)),      # Intercept
        x1,                    # x1 main effect
        x2,                    # x2 main effect
        x3,                    # x3 main effect
        x1 * x2,               # x1*x2 interaction
        x1 * x3,               # x1*x3 interaction
        x2 * x3,               # x2*x3 interaction
        x1 * x2 * x3           # x1*x2*x3 interaction
    ]

    # Use least squares to find the coefficients for the full model
    coeffs, _, _, _ = np.linalg.lstsq(X_full, y, rcond=None)

    # Define the names for each term in the model corresponding to the coefficients
    term_names = ["", "x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"]

    # Build the equation string for the parsimonious model
    # The intercept is always included in the model
    b0_rounded = round(coeffs[0], 1)
    equation_str = f"y = {b0_rounded}"

    # Iterate through the other terms and add them to the model if their effect is significant.
    # Based on the principle of parsimony, we define "significant" as having a 
    # coefficient magnitude greater than 1.0, which clearly separates the large from small effects.
    for i in range(1, len(coeffs)):
        coeff = coeffs[i]
        if abs(coeff) > 1.0:
            coeff_rounded = round(coeff, 1)
            term_name = term_names[i]
            
            # Determine the sign for the equation
            sign = "+" if coeff_rounded >= 0 else "-"
            value = abs(coeff_rounded)
            
            # Append ' sign value*term_name' to the equation string
            equation_str += f" {sign} {value}*{term_name}"
    
    print(equation_str)

create_parsimonious_model()