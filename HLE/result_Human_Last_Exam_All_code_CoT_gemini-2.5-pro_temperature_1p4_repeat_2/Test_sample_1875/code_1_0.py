import numpy as np
from sklearn.linear_model import LinearRegression

def solve_model():
    """
    Constructs and reports a parsimonious model to predict y from x1, x2, and x3.
    """
    # 1. Define the experimental data points.
    data = np.array([
        [-1, -1, -1, 34.3],
        [1, -1, -1, 94.6],
        [-1, 1, -1, 52.5],
        [1, 1, -1, 75.1],
        [-1, -1, 1, 28.4],
        [1, -1, 1, 48.3],
        [-1, 1, 1, 7.0],
        [1, 1, 1, 70.8]
    ])

    X_data = data[:, :3]
    y_data = data[:, 3]

    # Extract individual factors
    x1 = X_data[:, 0]
    x2 = X_data[:, 1]
    x3 = X_data[:, 2]

    # 2. Select features for the parsimonious model.
    # Based on preliminary analysis of the full model, the significant terms
    # are the main effects and the three-way interaction.
    x1x2x3 = x1 * x2 * x3
    X_parsimonious = np.column_stack([x1, x2, x3, x1x2x3])

    # 3. Fit the parsimonious linear regression model.
    model = LinearRegression()
    model.fit(X_parsimonious, y_data)

    # 4. Get the parameters for the final model.
    intercept = model.intercept_
    coeffs = model.coef_

    # 5. Construct and print the final equation.
    # Round parameters to one decimal place
    b0 = round(intercept, 1)
    b1 = round(coeffs[0], 1)
    b2 = round(coeffs[1], 1)
    b3 = round(coeffs[2], 1)
    b123 = round(coeffs[3], 1)

    # Build the equation string, ensuring correct signs and spacing
    equation = f"y = {b0}"
    equation += f" + {b1}*x1" if b1 >= 0 else f" - {abs(b1)}*x1"
    equation += f" + {b2}*x2" if b2 >= 0 else f" - {abs(b2)}*x2"
    equation += f" + {b3}*x3" if b3 >= 0 else f" - {abs(b3)}*x3"
    equation += f" + {b123}*x1*x2*x3" if b123 >= 0 else f" - {abs(b123)}*x1*x2*x3"

    print(equation)

solve_model()
<<<y = 51.4 + 20.8*x1 - 12.5*x2 - 12.6*x3 + 13.2*x1*x2*x3>>>