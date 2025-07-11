import numpy as np

def calculate_regression_coefficients():
    """
    Calculates multiple linear regression coefficients for the given epigenetic data.
    """
    # Step 1: Define the data
    # Independent variables (predictors)
    methylation = np.array([10, 15, 20, 25, 30, 35, 40])
    histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

    # Dependent variable (response)
    tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

    # Step 2: Prepare the data for regression
    # Create the feature matrix X with a column for the intercept
    # The intercept column (of ones) allows the model to calculate β0
    X = np.vstack([np.ones(len(methylation)), methylation, histone_h3k9]).T
    y = tsg_expression

    # Step 3 & 4: Calculate coefficients using numpy's least squares solver
    # This method is robust and handles the perfect multicollinearity in the data
    try:
        coefficients, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
        beta_0, beta_1, beta_2 = coefficients
    except np.linalg.LinAlgError as e:
        print(f"An error occurred during regression calculation: {e}")
        return

    # Step 5: Display the results
    print("Multiple Linear Regression Analysis Results:")
    print(f"Intercept (β0): {beta_0}")
    print(f"Coefficient for Methylation (β1): {beta_1}")
    print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

    print("\nThe regression equation is:")
    # Using '+' for negative numbers to ensure correct sign formatting
    print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9 Trimethylation)")


if __name__ == "__main__":
    calculate_regression_coefficients()

<<<(-100.0, 0.0, 2.0)>>>