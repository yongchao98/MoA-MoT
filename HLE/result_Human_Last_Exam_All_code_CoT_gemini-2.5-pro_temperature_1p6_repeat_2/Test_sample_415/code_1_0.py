import numpy as np
import statsmodels.api as sm

# --- 1. Define the Data ---
# Independent variables (Predictors)
# X1: Methylation (%), X2: Histone H3K9 Trimethylation
X = np.array([
    [10, 300],
    [15, 275],
    [20, 250],
    [25, 225],
    [30, 200],
    [35, 175],
    [40, 150]
])

# Dependent variable (Outcome)
# Y: TSG Expression (RNA-Seq Reads)
y = np.array([500, 450, 400, 350, 300, 250, 200])

# --- 2. Prepare Data and Fit the Model ---
# Add a constant (column of ones) to the predictor matrix to calculate the intercept
X = sm.add_constant(X)

# Create and fit the Ordinary Least Squares (OLS) model
model = sm.OLS(y, X)
results = model.fit()

# --- 3. Extract and Print Coefficients ---
# The results.params attribute holds the coefficients: [β0, β1, β2]
beta_0, beta_1, beta_2 = results.params

print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nFinal Regression Equation:")
# We use .2f to format the numbers to two decimal places for clarity
print(f"TSG Expression = {beta_0:.2f} + ({beta_1:.2f} * Methylation) + ({beta_2:.2f} * Histone H3K9 Trimethylation)")
