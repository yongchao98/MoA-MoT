import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Define the data from the table
# Independent variables (features)
methylation = np.array([10, 15, 20, 25, 30, 35, 40]).reshape(-1, 1)
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150]).reshape(-1, 1)

# Dependent variable (target)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# 2. Combine independent variables into a single feature matrix X
# The shape should be (n_samples, n_features) -> (7, 2)
X = np.hstack((methylation, histone_h3k9))
y = tsg_expression

# 3. Create and fit the multiple linear regression model
model = LinearRegression()
model.fit(X, y)

# 4. Extract the coefficients
beta_0 = model.intercept_
beta_1, beta_2 = model.coef_

# 5. Print the results
print("The multiple linear regression equation is:")
# The f-string formats the numbers to 4 decimal places for clarity
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone_H3K9_Trimethylation)")
print("\nCalculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0:.4f}")
print(f"Coefficient for Methylation (β1): {beta_1:.4f}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2:.4f}")
