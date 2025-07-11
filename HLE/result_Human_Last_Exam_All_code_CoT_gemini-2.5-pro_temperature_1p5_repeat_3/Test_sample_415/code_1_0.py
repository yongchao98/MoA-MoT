import numpy as np
from sklearn.linear_model import LinearRegression

# --- Data from the user's table ---
# Independent variable 1: Methylation (%)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])

# Independent variable 2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone_trimethylation = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# --- Multiple Linear Regression ---
# Combine independent variables into a single feature matrix X.
# Each row represents a sample, and each column represents a feature.
X = np.column_stack((methylation, histone_trimethylation))
y = tsg_expression

# Create a linear regression model instance
model = LinearRegression()

# Fit the model to the data
model.fit(X, y)

# --- Extract Coefficients ---
# b0 is the intercept of the regression plane
b0 = model.intercept_
# b1 and b2 are the slope coefficients for our independent variables
b1, b2 = model.coef_

# --- Print the Results ---
print("Multiple Linear Regression Results:")
print(f"β0 (Intercept): {b0}")
print(f"β1 (Coefficient for Methylation): {b1}")
print(f"β2 (Coefficient for Histone H3K9 Trimethylation): {b2}")
print("\n---")
print("The final regression equation is:")
print(f"TSG Expression = {b0:.2f} + ({b1:.2f} * Methylation) + ({b2:.2f} * Histone_H3K9_Trimethylation)")