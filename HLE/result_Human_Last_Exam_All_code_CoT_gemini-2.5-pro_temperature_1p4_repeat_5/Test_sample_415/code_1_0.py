import numpy as np
from sklearn.linear_model import LinearRegression

def calculate_regression_coefficients():
    """
    Calculates the multiple linear regression coefficients for the given dataset.
    """
    # 1. Define the data from the user's table
    methylation = np.array([10, 15, 20, 25, 30, 35, 40])
    histone = np.array([300, 275, 250, 225, 200, 175, 150])
    tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

    # 2. Prepare the data for the model
    # The independent variables (predictors) are combined into a single matrix X
    X = np.column_stack((methylation, histone))
    # The dependent variable (response) is y
    y = tsg_expression

    # 3. Create and fit the multiple linear regression model
    model = LinearRegression()
    model.fit(X, y)

    # 4. Extract the regression coefficients
    # β0 is the intercept
    beta_0 = model.intercept_
    # β1 and β2 are the slope coefficients for methylation and histone, respectively
    beta_1, beta_2 = model.coef_

    # 5. Print the results
    print("Multiple Linear Regression Results:")
    print(f"β0 (Intercept): {beta_0:.4f}")
    print(f"β1 (Coefficient for Methylation): {beta_1:.4f}")
    print(f"β2 (Coefficient for Histone H3K9 Trimethylation): {beta_2:.4f}")

    print("\nFinal Regression Equation:")
    # Print the equation with each coefficient value
    print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9 Trimethylation)")

# Execute the function to get the results
calculate_regression_coefficients()