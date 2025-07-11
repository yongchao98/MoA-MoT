import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data from the problem ---
# Independent variable 1: Methylation (%)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])

# Independent variable 2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Multiple Linear Regression ---
# Combine the independent variables into a single feature matrix (X)
# Each row represents a sample, and each column represents a feature.
X = np.column_stack((methylation, histone_h3k9))
y = tsg_expression

# Create a Linear Regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Extract and Print Coefficients ---
# Get the intercept (β0)
beta_0 = model.intercept_

# Get the coefficients for the independent variables (β1, β2)
beta_1, beta_2 = model.coef_

# Print the individual coefficients
print(f"Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")
print("\n" + "="*50 + "\n")

# Print the final regression equation
# Note: The coefficients are rounded for display purposes in the equation
print("Final Regression Equation:")
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9 Trimethylation)")