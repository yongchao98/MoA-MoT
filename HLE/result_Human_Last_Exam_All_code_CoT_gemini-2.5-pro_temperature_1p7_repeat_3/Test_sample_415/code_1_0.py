import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data ---
# Independent variables (features)
# Column 1: Methylation (%) at TSG Promoter
# Column 2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
X = np.array([
    [10, 300],
    [15, 275],
    [20, 250],
    [25, 225],
    [30, 200],
    [35, 175],
    [40, 150]
])

# Dependent variable
# TSG Expression (RNA-Seq Reads)
y = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Calculation ---
# Create a linear regression model object
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# Extract the coefficients
# b0 is the intercept
# b1 is the coefficient for the first feature (Methylation)
# b2 is the coefficient for the second feature (Histone H3K9 Trimethylation)
b0 = model.intercept_
b1, b2 = model.coef_

# --- Output ---
print("Multiple Linear Regression Coefficients:")
print(f"Intercept (β0): {b0}")
print(f"Coefficient for Methylation (β1): {b1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {b2}")
print("\nFinal Regression Equation:")
print(f"TSG Expression = {b0:.4f} + ({b1:.4f} * Methylation) + ({b2:.4f} * Histone_H3K9_Trimethylation)")

final_answer = f"β0 = {b0}, β1 = {b1}, β2 = {b2}"
# The <<<...>>> format is for machine readability, the formatted text above is for human readability.
# Do not copy the line below. It is used for auto-extraction.
print(f'<<<{final_answer}>>>')