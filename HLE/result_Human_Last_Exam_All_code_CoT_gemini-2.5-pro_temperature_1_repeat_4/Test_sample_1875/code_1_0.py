import numpy as np

def solve_model_and_predict():
    """
    Constructs and presents a parsimonious model for the given data.
    """
    # 1. Define the data
    # The data consists of 8 points from a 2^3 factorial design.
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])

    # 2. Formulate the full model and estimate parameters
    # The design matrix X includes columns for the intercept, main effects,
    # and all interaction terms.
    X_full = np.column_stack([
        np.ones(8),      # Intercept
        x1,              # x1
        x2,              # x2
        x3,              # x3
        x1 * x2,         # x1*x2 interaction
        x1 * x3,         # x1*x3 interaction
        x2 * x3,         # x2*x3 interaction
        x1 * x2 * x3     # x1*x2*x3 interaction
    ])
    
    term_labels = [
        '', ' * x1', ' * x2', ' * x3', 
        ' * x1 * x2', ' * x1 * x3', ' * x2 * x3', ' * x1 * x2 * x3'
    ]

    # Use least squares to solve for the coefficients b in y = Xb
    coeffs, _, _, _ = np.linalg.lstsq(X_full, y, rcond=None)

    # 3. Identify the least significant term for a parsimonious model
    # The term with the coefficient of the smallest magnitude is the least significant.
    # We exclude the intercept (b0) from this search.
    magnitudes = np.abs(coeffs[1:])
    min_mag_index = np.argmin(magnitudes) + 1 # Add 1 to account for intercept
    
    # In this case, the coefficient for x1*x2 (coeffs[4] = 0.775) is the smallest.
    # We will remove this term to create the parsimonious model.
    # Due to orthogonality, other coefficients remain unchanged.
    
    # 4. Construct and print the final model equation
    # Start the equation with the intercept term.
    equation = f"y = {coeffs[0]:.1f}"
    
    # Append the other significant terms to the equation string.
    for i in range(1, len(coeffs)):
        # Skip the term being removed.
        if i == min_mag_index:
            continue
            
        coeff_val = coeffs[i]
        term_str = term_labels[i]
        
        # Format the string with the correct sign and value.
        sign = '+' if coeff_val >= 0 else '-'
        equation += f" {sign} {abs(coeff_val):.1f}{term_str}"
        
    print("The parsimonious model is:")
    print(equation)
    
    # Also print the parameter estimates clearly
    print("\nParameter estimates (rounded to one decimal place):")
    param_names = ['b0 (Intercept)', 'b1 (for x1)', 'b2 (for x2)', 'b3 (for x3)', 
                   'b12 (for x1*x2)', 'b13 (for x1*x3)', 'b23 (for x2*x3)', 'b123 (for x1*x2*x3)']
    
    final_params = []
    for i in range(len(coeffs)):
        if i == min_mag_index:
            line = f"{param_names[i]}: (Removed from model)"
        else:
            line = f"{param_names[i]}: {coeffs[i]:.1f}"
        final_params.append(line)
        print(line)

solve_model_and_predict()
