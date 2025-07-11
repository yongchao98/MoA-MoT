import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Load the Data
# Independent variables from the table
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable from the table
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# 2. Define Variables for the model
# The independent variables (features) are combined into a single 2D array (matrix)
# Each row represents a sample, and each column represents a feature.
X = np.array([methylation, histone_h3k9]).T

# The dependent variable (target)
y = tsg_expression

# 3. Perform Regression & 4. Fit the Model
# Create a Linear Regression model instance
model = LinearRegression()

# Fit the model to the data to calculate the coefficients
model.fit(X, y)

# 5. Extract Coefficients
# The intercept (β0)
beta_0 = model.intercept_
# The slope coefficients for methylation (β1) and histone modification (β2)
beta_1, beta_2 = model.coef_

# 6. Display Results
print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nFinal Regression Equation:")
# We use .4f to format the numbers to 4 decimal places for readability in the equation string
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f}) * Methylation + ({beta_2:.4f}) * Histone_H3K9_Trimethylation")
