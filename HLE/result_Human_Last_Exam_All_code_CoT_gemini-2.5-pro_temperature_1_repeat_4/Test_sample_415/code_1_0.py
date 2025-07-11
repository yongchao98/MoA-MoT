import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data from the user's table ---
# Independent variables (predictors)
# X1: Methylation (%) at TSG Promoter
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
# X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (response)
# Y: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Multiple Linear Regression ---
# We need to combine our independent variables into a single feature matrix.
# The shape will be (n_samples, n_features), which is (7, 2) in this case.
X = np.stack((methylation, histone_h3k9), axis=1)
y = tsg_expression

# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Extract Coefficients ---
# β0 is the intercept
beta_0 = model.intercept_
# β1 and β2 are the slope coefficients for our two features
beta_1, beta_2 = model.coef_

# --- Print the Results ---
print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nFinal Regression Equation:")
# The equation is: Y = β0 + β1*X1 + β2*X2
print(f"TSG Expression = {beta_0:.2f} + ({beta_1:.2f} * Methylation) + ({beta_2:.2f} * Histone H3K9 Trimethylation)")

print(f"\n<<<β0 = {beta_0}, β1 = {beta_1}, β2 = {beta_2}>>>")