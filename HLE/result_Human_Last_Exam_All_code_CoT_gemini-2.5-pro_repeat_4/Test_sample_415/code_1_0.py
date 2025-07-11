import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Preparation ---
# Independent variables (features)
# X1: Methylation (%)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
# X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (target)
# Y: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# For scikit-learn, combine independent variables into a single feature matrix
X = np.column_stack((methylation, histone_h3k9))
y = tsg_expression

# --- Model Training ---
# Create an instance of the Linear Regression model
model = LinearRegression()

# Fit the model to the data to find the best-fit line
model.fit(X, y)

# --- Coefficient Extraction ---
# Get the intercept (β0) of the regression line
beta_0 = model.intercept_

# Get the slope coefficients (β1 for methylation, β2 for histone H3K9)
beta_1, beta_2 = model.coef_

# --- Output Results ---
print("The calculated regression coefficients are:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nThe final multiple linear regression equation is:")
# The format below prints each number in the final equation.
# β2 is extremely close to zero, which is common with collinear data.
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9)")
