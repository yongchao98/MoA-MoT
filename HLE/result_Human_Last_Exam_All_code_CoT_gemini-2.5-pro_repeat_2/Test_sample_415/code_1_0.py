import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Load Data
# Independent variable X1: Methylation (%)
methylation = np.array([10, 15, 20, 25, 30, 35, 40]).reshape(-1, 1)

# Independent variable X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone_h3k9 = np.array([300, 275, 250, 225, 200, 175, 150]).reshape(-1, 1)

# Dependent variable Y: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Combine independent variables into a feature matrix
X = np.hstack((methylation, histone_h3k9))
y = tsg_expression

# 2. Set up and 3. Fit the Model
model = LinearRegression()
model.fit(X, y)

# 4. Extract and Display Coefficients
beta_0 = model.intercept_
beta_1, beta_2 = model.coef_

print("Multiple Linear Regression Results:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")
print("\n---")
print("The final regression equation is:")
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9 Trimethylation)")

# The data in this problem exhibits perfect multicollinearity, meaning the independent variables
# are perfectly correlated. While this can make individual coefficients hard to interpret in a
# biological context as their effects are confounded, the model can still be used for prediction.
# The scikit-learn library provides one of the possible solutions for the coefficients.

# Final answer in the requested format
final_coeffs = f"{beta_0},{beta_1},{beta_2}"
# print(f"\n<<<{final_coeffs}>>>")