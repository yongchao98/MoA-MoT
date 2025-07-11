import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data from the user's table ---

# Independent variables (predictors)
# X1: Methylation (%) at TSG Promoter
# X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
X_data = np.array([
    [10, 300],
    [15, 275],
    [20, 250],
    [25, 225],
    [30, 200],
    [35, 175],
    [40, 150]
])

# Dependent variable (response)
# Y: TSG Expression (RNA-Seq Reads)
y_data = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Multiple Linear Regression ---

# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X_data, y_data)

# Extract the regression coefficients
beta_0 = model.intercept_  # Intercept (β0)
beta_1, beta_2 = model.coef_  # Slope coefficients (β1 for Methylation, β2 for Histone)

# --- Output the Results ---

print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")
print("\nFinal Regression Equation:")
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9 Trimethylation)")

# The final answer for the platform
final_answer = f"{beta_0}, {beta_1}, {beta_2}"
# print(f"\n<<< {final_answer} >>>") # This is for generating the final answer string