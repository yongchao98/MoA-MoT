import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Preparation ---
# Independent variables (predictors) based on the provided table.
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

# Dependent variable (outcome)
# y: TSG Expression (RNA-Seq Reads)
y = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Model Training ---
# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data to calculate the coefficients
model.fit(X, y)

# --- Extract Coefficients ---
# Intercept (β0)
beta_0 = model.intercept_

# Slope coefficients (β1 for Methylation, β2 for Histone)
beta_1, beta_2 = model.coef_

# --- Output Results ---
# The data exhibits perfect multicollinearity, meaning the independent variables
# are perfectly correlated. The scikit-learn library provides one of the possible
# valid solutions for the coefficients.

print("The calculated multiple linear regression equation is:")
# The format string handles positive/negative signs correctly for the equation
# Note: The coefficient for Histone is effectively zero due to multicollinearity.
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9 Trimethylation)")
print("\n" + "-"*25)
print("The regression coefficients are:")
print(f"β0 (Intercept): {beta_0}")
print(f"β1 (Coefficient for Methylation): {beta_1}")
print(f"β2 (Coefficient for Histone H3K9 Trimethylation): {beta_2}")