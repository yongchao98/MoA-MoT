import numpy as np
from sklearn.linear_model import LinearRegression

# --- 1. Data Preparation ---
# Independent variables (predictors)
# X1: Methylation (%)
# X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
X = np.array([
    [10, 300],
    [15, 275],
    [20, 250],
    [25, 225],
    [30, 200],
    [35, 175],
    [40, 150]
])

# Dependent variable (outcome)
# y: TSG Expression (RNA-Seq Reads)
y = np.array([500, 450, 400, 350, 300, 250, 200])

# --- 2. Model Selection & 3. Fitting ---
# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- 4. Coefficient Extraction ---
# Get the intercept (beta_0)
beta_0 = model.intercept_

# Get the slope coefficients (beta_1 for Methylation, beta_2 for Histone)
beta_1, beta_2 = model.coef_

# --- Output ---
# Print the final regression equation with the calculated coefficients
print("The calculated multiple linear regression equation is:")
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f}) * Methylation + {beta_2:.4f} * Histone_H3K9_Trimethylation")

print("\nWhere the coefficients are:")
print(f"β0 (Intercept): {beta_0:.4f}")
print(f"β1 (Methylation): {beta_1:.4f}")
print(f"β2 (Histone H3K9 Trimethylation): {beta_2:.4f}")