import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Preparation ---
# Independent variables (features): 
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

# Dependent variable (target): 
# Y: TSG Expression (RNA-Seq Reads)
y = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Model Fitting ---
# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Extract Coefficients ---
# b0 is the intercept of the regression plane
beta_0 = model.intercept_

# b1 and b2 are the coefficients for the independent variables
beta_1, beta_2 = model.coef_

# --- Display Results ---
print("Multiple Linear Regression Results:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nThe final regression equation is:")
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f}) * Methylation + ({beta_2:.4f}) * Histone_H3K9_Trimethylation")

# The format requested by the prompt is to return the values in the final answer
final_answer = f"β0={beta_0}, β1={beta_1}, β2={beta_2}"
print(f"\n<<<{final_answer}>>>")