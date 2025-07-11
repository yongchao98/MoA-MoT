import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Preparation ---
# Independent variables (predictors)
# X1: Methylation (%)
# X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
X = np.array([
    [10, 300],  # Sample 1
    [15, 275],  # Sample 2
    [20, 250],  # Sample 3
    [25, 225],  # Sample 4
    [30, 200],  # Sample 5
    [35, 175],  # Sample 6
    [40, 150]   # Sample 7
])

# Dependent variable (response)
# y: TSG Expression (RNA-Seq Reads)
y = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Model Training ---
# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Results ---
# Extract the coefficients from the trained model
# b0 is the intercept of the regression line
b0 = model.intercept_
# b1 and b2 are the slope coefficients for methylation and histone trimethylation, respectively
b1, b2 = model.coef_

# --- Output ---
print("The calculated multiple linear regression equation is:")
print(f"TSG Expression = {b0:.4f} + ({b1:.4f} * Methylation) + ({b2:.4f} * Histone H3K9 Trimethylation)")
print("\n---")
print("Regression Coefficients:")
print(f"Intercept (β0): {b0:.4f}")
print(f"Methylation Coefficient (β1): {b1:.4f}")
print(f"Histone H3K9 Trimethylation Coefficient (β2): {b2:.4f}")
