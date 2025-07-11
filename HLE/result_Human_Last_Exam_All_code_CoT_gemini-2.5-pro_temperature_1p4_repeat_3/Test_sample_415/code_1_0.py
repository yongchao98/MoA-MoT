import numpy as np
from sklearn.linear_model import LinearRegression

def calculate_epigenetic_regression():
    """
    Calculates the multiple linear regression for the provided epigenetic data.
    """
    # Independent variables (features)
    # Methylation (%) at TSG Promoter
    methylation = np.array([10, 15, 20, 25, 30, 35, 40])
    # Histone H3K9 Trimethylation (ChIP-Seq Reads)
    histone = np.array([300, 275, 250, 225, 200, 175, 150])

    # Dependent variable (target)
    # TSG Expression (RNA-Seq Reads)
    tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

    # We need to combine the independent variables into a single feature matrix.
    # Each row is a sample, each column is a feature.
    X = np.column_stack((methylation, histone))
    y = tsg_expression

    # Create and fit the multiple linear regression model
    model = LinearRegression()
    model.fit(X, y)

    # Extract the coefficients
    beta_0 = model.intercept_
    beta_1, beta_2 = model.coef_

    # Print the individual coefficients
    print("Multiple Linear Regression Coefficients:")
    print(f"Intercept (β0): {beta_0:.4f}")
    print(f"Coefficient for Methylation (β1): {beta_1:.4f}")
    print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2:.4f}")

    # Print the final formatted equation
    print("\nResulting Equation:")
    # We use round() for a cleaner display and check the sign for correct formatting
    b0_disp = round(beta_0, 2)
    b1_disp = round(beta_1, 2)
    b2_disp = round(beta_2, 2)
    
    sign_beta1 = "-" if b1_disp < 0 else "+"
    sign_beta2 = "-" if b2_disp < 0 else "+"

    print(f"TSG Expression = {b0_disp} {sign_beta1} {abs(b1_disp)} * Methylation {sign_beta2} {abs(b2_disp)} * Histone_H3K9_Trimethylation")

if __name__ == '__main__':
    calculate_epigenetic_regression()