import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Define the data from the table
# Independent variables (features)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])  # X1
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150]) # X2

# Dependent variable (target)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200]) # Y

# 2. Set up the regression model
# Scikit-learn requires the independent variables to be in a 2D array format
X = np.column_stack((methylation, histone_h3k9))
y = tsg_expression

# 3. Create a Linear Regression model object and fit it to the data
model = LinearRegression()
model.fit(X, y)

# 4. Extract and display the coefficients
beta_0 = model.intercept_
beta_1, beta_2 = model.coef_

print("Multiple Linear Regression Results:")
print(f"Intercept (β0): {beta_0:.4f}")
print(f"Coefficient for Methylation (β1): {beta_1:.4f}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2:.4f}")

print("\nThe final regression equation is:")
print(f"TSG Expression = {beta_0:.2f} + ({beta_1:.2f} * Methylation) + ({beta_2:.2f} * Histone_H3K9_Trimethylation)")
