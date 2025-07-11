import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Preparation ---
# Independent variables (predictors)
# X1: Methylation (%)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
# X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (outcome)
# Y: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Combine independent variables into a single feature matrix for scikit-learn
# The .T transposes the matrix to the correct shape (samples, features)
X = np.vstack((methylation, histone_h3k9)).T
y = tsg_expression

# --- Model Training ---
# Create a Linear Regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Results ---
# Get the regression coefficients
# β0: Intercept
beta_0 = model.intercept_
# β1, β2: Slope coefficients for methylation and histone modification
beta_1, beta_2 = model.coef_

# Print the results in the final equation format
print("Multiple Linear Regression Results:")
print("The calculated regression equation is:")
print(f"TSG Expression = {beta_0:.1f} + ({beta_1:.1f}) * Methylation + ({beta_2:.1f}) * Histone_H3K9_Trimethylation")
print("\nIndividual Coefficients:")
print(f"Intercept (β0): {beta_0:.1f}")
print(f"Coefficient for Methylation (β1): {beta_1:.1f}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2:.1f}")
