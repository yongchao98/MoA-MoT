import numpy as np
from sklearn.linear_model import LinearRegression

# Step 1: Define the data based on the provided table
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Step 2: Construct the independent variable matrix (X)
# We combine the two independent variables into a single matrix.
X = np.vstack((methylation, histone_h3k9)).T

# We add a constant column of ones to X for the intercept calculation.
X_with_const = np.c_[np.ones(X.shape[0]), X]

# Step 3: Set up the dependent variable vector (y)
y = tsg_expression

# Step 4: Create and fit the multiple linear regression model
# We set fit_intercept=False because we've already added the constant column.
model = LinearRegression(fit_intercept=False)
model.fit(X_with_const, y)

# Step 6: Extract the coefficients
# The coefficients will be in the order they appear in X_with_const: [intercept, methylation, histone]
beta_0, beta_1, beta_2 = model.coef_

# Step 7: Print the results
print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nFinal Regression Equation:")
# We use .2f to format the numbers to two decimal places for clarity.
print(f"TSG Expression = {beta_0:.2f} + ({beta_1:.2f} * Methylation) + ({beta_2:.2f} * Histone_H3K9_Trimethylation)")
