import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Preparation ---
# Independent variables (features) from the table
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (target) from the table
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Combine independent variables into a single feature matrix.
# The shape will be (7 samples, 2 features).
X = np.column_stack((methylation, histone_h3k9))
y = tsg_expression

# --- Model Training ---
# Create and fit the multiple linear regression model
model = LinearRegression()
model.fit(X, y)

# --- Results ---
# Extract the coefficients from the fitted model
beta_0 = model.intercept_
beta_1, beta_2 = model.coef_

# Print the calculated coefficients
print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0:.4f}")
print(f"Coefficient for Methylation (β1): {beta_1:.4f}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2:.4f}")
print("\n---")

# Print the final regression equation as requested
print("Final Regression Equation:")
print(f"TSG Expression = {beta_0:.2f} + ({beta_1:.2f}) * Methylation + ({beta_2:.2f}) * Histone_H3K9")

# Note: The data exhibits perfect multicollinearity, meaning one independent variable
# is a perfect linear function of the other. Different statistical packages might
# attribute the effect differently. The standard sklearn result is provided here,
# where the effect is primarily attributed to one variable.