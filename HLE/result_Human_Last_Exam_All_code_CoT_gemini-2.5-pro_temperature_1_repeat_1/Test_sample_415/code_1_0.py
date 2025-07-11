import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data from the user ---
# Independent variables (predictors)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (outcome)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Prepare data for scikit-learn ---
# Combine the independent variables into a single feature matrix X
# Each row is a sample, each column is a feature
X = np.vstack((methylation, histone_h3k9)).T
y = tsg_expression

# --- Perform Multiple Linear Regression ---
# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Extract the coefficients ---
# beta_0 is the intercept of the regression plane
beta_0 = model.intercept_

# beta_1 and beta_2 are the coefficients for methylation and histone H3K9, respectively
beta_1, beta_2 = model.coef_

# --- Display the results ---
print("Multiple Linear Regression Calculation:")
print("-" * 40)
print(f"The calculated regression equation is:")
print(f"TSG Expression = {beta_0:.2f} + ({beta_1:.2f} * Methylation) + ({beta_2:.2f} * Histone H3K9)")
print("-" * 40)
print("The regression coefficients are:")
print(f"Intercept (β0): {beta_0:.4f}")
print(f"Coefficient for Methylation (β1): {beta_1:.4f}")
print(f"Coefficient for Histone H3K9 (β2): {beta_2:.4f}")

# Final answer format
final_answer = f"β0 = {beta_0:.4f}, β1 = {beta_1:.4f}, β2 = {beta_2:.4f}"
print(f"\n<<<{final_answer}>>>")