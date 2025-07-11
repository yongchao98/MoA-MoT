import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data Preparation ---
# Independent variables (features) from your table
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_mod = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (target) from your table
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# We need to combine our independent variables into a single 2D array (feature matrix)
# where each row is a sample and each column is a feature.
X = np.column_stack((methylation, histone_mod))
Y = tsg_expression

# --- Model Calculation ---
# Create a Linear Regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, Y)

# --- Results ---
# Extract the coefficients from the fitted model
beta_0 = model.intercept_
coefficients = model.coef_
beta_1 = coefficients[0]
beta_2 = coefficients[1]

# Print the individual coefficients
print("Calculated Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

# Print the final regression equation
# Determine the signs for a clean equation string
sign_beta_1 = "-" if beta_1 < 0 else "+"
sign_beta_2 = "-" if beta_2 < 0 else "+"

print("\nFinal Regression Equation:")
# The following line prints the full equation with the calculated coefficients
print(f"TSG Expression = {beta_0:.4f} {sign_beta_1} {abs(beta_1):.4f} * Methylation {sign_beta_2} {abs(beta_2):.4f} * Histone_H3K9_Trimethylation")

# Output the final answer in the specified format
final_answer = (beta_0, beta_1, beta_2)
# print(f"\n<<<{final_answer}>>>") # This is for generating the final response below.