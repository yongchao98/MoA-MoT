import numpy as np

def solve_parsimonious_model():
    """
    Constructs and prints a parsimonious model for the given dataset.
    """
    # The available data points
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

    x_data = data[:, :3]
    y_data = data[:, 3]

    # 1. Construct the design matrix for the full model including all interactions
    X_full_design = np.c_[
        np.ones(x_data.shape[0]),      # Intercept
        x_data[:, 0],                  # x1
        x_data[:, 1],                  # x2
        x_data[:, 2],                  # x3
        x_data[:, 0] * x_data[:, 1],    # x1*x2
        x_data[:, 0] * x_data[:, 2],    # x1*x3
        x_data[:, 1] * x_data[:, 2],    # x2*x3
        x_data[:, 0] * x_data[:, 1] * x_data[:, 2] # x1*x2*x3
    ]
    
    term_names = ["", "x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"]

    # 2. Solve for the coefficients using least squares
    coeffs, _, _, _ = np.linalg.lstsq(X_full_design, y_data, rcond=None)

    # 3. Identify significant terms by magnitude.
    # Coefficients for x2 (index 2) and x1*x2 (index 4) are very small.
    # We will exclude them to create a parsimonious model.
    
    # 4. Define the terms and coefficients for the parsimonious model
    parsimonious_terms = {
        "": coeffs[0],
        "x1": coeffs[1],
        "x3": coeffs[3],
        "x1*x3": coeffs[5],
        "x2*x3": coeffs[6],
        "x1*x2*x3": coeffs[7]
    }
    
    # Define the order for a readable equation
    term_order = ["", "x1", "x3", "x1*x3", "x2*x3", "x1*x2*x3"]

    # 5. Construct the final equation string
    equation_parts = []
    
    # Handle the intercept first
    intercept_val = round(parsimonious_terms[""], 1)
    equation_parts.append(f"{intercept_val:.1f}")

    # Handle the other terms
    for term in term_order[1:]:
        coeff = parsimonious_terms[term]
        coeff_rounded = round(coeff, 1)
        
        # Determine the sign and format the term
        sign = "-" if coeff_rounded < 0 else "+"
        
        # Format the variable part of the term
        var_part = f"*{term}"
        
        equation_parts.append(f"{sign} {abs(coeff_rounded):.1f}{var_part}")

    final_equation = "y = " + " ".join(equation_parts)
    
    print(final_equation)

solve_parsimonious_model()
<<<y = 51.4 + 20.8*x1 - 12.6*x3 - 17.6*x1*x3 + 10.6*x2*x3 + 10.2*x1*x2*x3>>>