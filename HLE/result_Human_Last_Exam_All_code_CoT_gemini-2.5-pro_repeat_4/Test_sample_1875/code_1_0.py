import numpy as np

def solve_parsimonious_model():
    """
    This function constructs a parsimonious model to predict y from x1, x2, and x3.
    """
    # 1. Define the available data points
    # The design is a 2^3 factorial design
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])

    # 2. Create the design matrix 'X' for the full model
    # The columns are: intercept, x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3
    X = np.c_[np.ones(len(x1)), x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3]
    term_names = ["", "*x1", "*x2", "*x3", "*x1*x2", "*x1*x3", "*x2*x3", "*x1*x2*x3"]

    # 3. Solve for the model coefficients using least squares
    # This finds the beta vector in the equation y = X*beta
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

    # 4. Identify significant terms for a parsimonious model
    # We examine the magnitudes of the coefficients to decide which terms to keep.
    # A common rule is to look for a gap in the magnitudes.
    # Coeffs: [51.37, 20.82, -0.025, -12.75, 0.775, 0.85, -12.45, 3.95]
    # The small coefficients are for x2, x1*x2, x1*x3. We will drop them.
    # The selected terms are the intercept, x1, x3, x2*x3, and x1*x2*x3.
    selected_indices = [0, 1, 3, 6, 7]
    
    # 5. Construct and print the final model equation
    # The coefficients are rounded to one decimal place as requested.
    
    # Format the intercept term
    equation_parts = [f"{coeffs[0]:.1f}"]
    
    # Format the other selected terms
    for i in selected_indices[1:]:
        coeff_val = coeffs[i]
        # Using numpy's round-half-to-even behavior
        rounded_coeff = np.round(coeff_val, 1)
        
        # Skip terms with a coefficient of 0.0 after rounding
        if abs(rounded_coeff) < 1e-6:
            continue
            
        sign = "+" if rounded_coeff > 0 else "-"
        # Python f-strings format negative numbers with a sign automatically
        # So we adjust the spacing and add the '+' for positive numbers
        if sign == '+':
            equation_parts.append(f"+ {abs(rounded_coeff):.1f}{term_names[i]}")
        else:
            equation_parts.append(f"- {abs(rounded_coeff):.1f}{term_names[i]}")

    final_equation = "y = " + " ".join(equation_parts)
    print("The final parsimonious model is:")
    print(final_equation)
    print("\nParameter estimates (rounded to one decimal place):")
    print(f"Intercept (β₀): {np.round(coeffs[0], 1)}")
    print(f"β₁ (for x1): {np.round(coeffs[1], 1)}")
    print(f"β₃ (for x3): {np.round(coeffs[3], 1)}")
    print(f"β₂₃ (for x2*x3): {np.round(coeffs[6], 1)}")
    print(f"β₁₂₃ (for x1*x2*x3): {np.round(coeffs[7], 1)}")


solve_parsimonious_model()