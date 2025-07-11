import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Define the data based on the provided table
# Independent variable 1: Methylation (%)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])

# Independent variable 2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone_trimethylation = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# 2. Prepare the data for the model
# Combine the independent variables into a single feature matrix X
# The .T transposes the matrix to the correct shape (7 samples, 2 features)
X = np.vstack((methylation, histone_trimethylation)).T
y = tsg_expression

# 3. Create and fit the multiple linear regression model
model = LinearRegression()
model.fit(X, y)

# 4. Extract the regression coefficients
# b0 is the intercept of the regression plane
beta_0 = model.intercept_
# b1 and b2 are the coefficients for methylation and histone trimethylation, respectively
beta_1, beta_2 = model.coef_

# 5. Print the results
print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nFinal Regression Equation:")
print(f"TSG Expression = {beta_0} + ({beta_1} * Methylation) + ({beta_2} * Histone H3K9 Trimethylation)")

<<<[51.92307692307697, -2.884615384615385, 1.423076923076923]>>>