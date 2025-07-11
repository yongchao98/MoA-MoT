import numpy as np
from sklearn.linear_model import LinearRegression

# The provided data for the independent variables (Methylation and Histone Trimethylation)
# are perfectly correlated. In such a case of perfect multicollinearity, there are infinite
# valid solutions for the coefficients. A standard method to find a unique and stable solution
# is to remove the redundant variable. We will model the expression based on methylation and
# set the coefficient for the redundant histone variable to zero.

# 1. Define the primary independent variable (X) and the dependent variable (y)
methylation = np.array([10, 15, 20, 25, 30, 35, 40]).reshape(-1, 1)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Note: The histone data is perfectly redundant given the methylation data.
# The relationship is: Histone = 350 - 5 * Methylation.

# 2. Create and fit a simple linear regression model using only methylation
model = LinearRegression()
model.fit(methylation, tsg_expression)

# 3. Extract the coefficients for the simplified model.
# The intercept is our β0, and the methylation coefficient is our β1.
beta_0 = model.intercept_
beta_1 = model.coef_[0]

# The coefficient for the redundant histone variable, β2, is 0 in this solution.
beta_2 = 0.0

# 4. Print the final coefficients and the full regression equation
print("The calculated regression coefficients are:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nFinal Regression Equation:")
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f}) * Methylation + ({beta_2:.4f}) * Histone_H3K9_Trimethylation")