import numpy as np
from sklearn.linear_model import LinearRegression

# Step 1: Organize the data from the table
# Independent variables (features)
methylation = np.array([10, 15, 20, 25, 30, 35, 40]).reshape(-1, 1)
histone = np.array([300, 275, 250, 225, 200, 175, 150]).reshape(-1, 1)

# Dependent variable (target)
expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Step 2: Prepare the data for the regression model
# Combine the independent variables into a single feature matrix X
X = np.concatenate([methylation, histone], axis=1)
y = expression

# Step 3: Create and fit the multiple linear regression model
model = LinearRegression()
model.fit(X, y)

# Step 4: Extract the coefficients
# b0 is the intercept
beta_0 = model.intercept_
# b1 and b2 are the slope coefficients for methylation and histone modification
beta_1, beta_2 = model.coef_

# Step 5: Display the results
print("Multiple Linear Regression Results:")
print("-" * 35)
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")
print("-" * 35)
print("\nFinal Regression Equation:")
# The format ensures the signs of the coefficients are displayed correctly.
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone)")

final_answer_string = f"β0 = {beta_0}, β1 = {beta_1}, β2 = {beta_2}"
print(f"\n<<<{final_answer_string}>>>")