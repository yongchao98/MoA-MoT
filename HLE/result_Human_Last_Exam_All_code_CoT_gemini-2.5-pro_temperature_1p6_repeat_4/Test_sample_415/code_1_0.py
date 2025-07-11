import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Define the data from the table
# Independent variables (predictors)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_trimethylation = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (outcome)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# 2. Prepare the data for the model
# The scikit-learn library expects the independent variables (X) as a 2D array.
X = np.vstack([methylation, histone_trimethylation]).T
y = tsg_expression

# 3. Create and fit the multiple linear regression model
model = LinearRegression()
model.fit(X, y)

# 4. Extract the coefficients from the fitted model
beta_0 = model.intercept_
beta_1, beta_2 = model.coef_

# 5. Print the results
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")
print("\nThe final regression equation is:")
print(f"TSG Expression = {beta_0:.1f} + ({beta_1:.1f} * Methylation) + ({beta_2:.1f} * Histone H3K9 Trimethylation)")

# Final answer in the specified format for automated reading.
# The values are rounded to 4 decimal places for precision in the final answer block.
final_answer = f"<<<{beta_0:.4f}, {beta_1:.4f}, {beta_2:.4f}>>>"