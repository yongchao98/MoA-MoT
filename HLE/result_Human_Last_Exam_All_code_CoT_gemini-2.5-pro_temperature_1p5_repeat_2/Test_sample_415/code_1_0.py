import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data from the user's table ---
# Independent variables (predictors)
methylation = np.array([10, 15, 20, 25, 30, 35, 40]).reshape(-1, 1)
histone_mod = np.array([300, 275, 250, 225, 200, 175, 150]).reshape(-1, 1)

# Dependent variable (outcome)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Prepare data for multiple regression ---
# Combine the two independent variables into a single feature matrix
X = np.hstack([methylation, histone_mod])
y = tsg_expression

# --- Perform Multiple Linear Regression ---
# Create a linear regression model object
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Extract and display the coefficients ---
# β0 (Intercept)
b0 = model.intercept_
# β1 and β2 (Slope coefficients for the two variables)
b1, b2 = model.coef_

print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {b0:.2f}")
print(f"Coefficient for Methylation (β1): {b1:.2f}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {b2:.2f}")

# Display the final regression equation
print("\nFinal Regression Equation:")
print(f"TSG Expression = {b0:.2f} + ({b1:.2f} * Methylation) + ({b2:.2f} * Histone_H3K9_Trimethylation)")
