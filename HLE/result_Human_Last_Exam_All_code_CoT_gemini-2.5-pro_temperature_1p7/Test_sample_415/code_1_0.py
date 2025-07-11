import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Setup ---
# The data is taken from the table provided in the problem.

# Independent variables (Predictors):
# X1: Methylation (%) at TSG Promoter
# X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
X = np.array([
    [10, 300],  # Sample 1
    [15, 275],  # Sample 2
    [20, 250],  # Sample 3
    [25, 225],  # Sample 4
    [30, 200],  # Sample 5
    [35, 175],  # Sample 6
    [40, 150]   # Sample 7
])

# Dependent variable (Response):
# Y: TSG Expression (RNA-Seq Reads)
y = np.array([500, 450, 400, 350, 300, 250, 200])


# --- Multiple Linear Regression ---
# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data to calculate the coefficients
model.fit(X, y)

# Extract the calculated coefficients
# beta_0 is the intercept of the regression plane
beta_0 = model.intercept_
# beta_1 and beta_2 are the slope coefficients for methylation and histone modification, respectively
beta_1, beta_2 = model.coef_


# --- Output the Results ---
# Print the results in a clear and understandable format
print("Multiple Linear Regression Results")
print("-" * 40)
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")
print("-" * 40)
print("Final Regression Equation:")
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone_H3K9_Trimethylation)")
