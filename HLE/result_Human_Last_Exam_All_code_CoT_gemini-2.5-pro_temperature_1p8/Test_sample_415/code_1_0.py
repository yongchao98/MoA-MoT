import numpy as np
from sklearn.linear_model import LinearRegression

def calculate_regression_coefficients():
    """
    Calculates and prints multiple linear regression coefficients for gene expression data.
    """
    # Independent variables (predictors)
    # X1: Methylation (%)
    # X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
    X = np.array([
        [10, 300],
        [15, 275],
        [20, 250],
        [25, 225],
        [30, 200],
        [35, 175],
        [40, 150]
    ])

    # Dependent variable (response)
    # y: TSG Expression (RNA-Seq Reads)
    y = np.array([500, 450, 400, 350, 300, 250, 200])

    # Create and fit the multiple linear regression model
    model = LinearRegression()
    model.fit(X, y)

    # Extract the coefficients
    beta_0 = model.intercept_
    beta_1, beta_2 = model.coef_

    # Print the final regression equation
    print("Combined Effect Model (Multiple Linear Regression):")
    # Using '+' for negative coefficients to clearly show the equation structure
    print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone_H3K9_Trimethylation)\n")

    # Print the individual coefficients
    print("Regression Coefficients:")
    print(f"Intercept (β0): {beta_0:.4f}")
    print(f"Coefficient for Methylation (β1): {beta_1:.4f}")
    print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2:.4f}")
    
    # Return a tuple for the final answer format
    return (beta_0, beta_1, beta_2)

# Run the calculation and get the coefficients for the final answer
final_coeffs = calculate_regression_coefficients()

# The final answer format is handled outside the printed output
# for submission systems, but let's make it explicit here for clarity.
final_answer_string = f"<<<{final_coeffs}>>>"
# In a real script, we would just use the values directly.
# For this task, the final answer will be formatted at the end of the response.