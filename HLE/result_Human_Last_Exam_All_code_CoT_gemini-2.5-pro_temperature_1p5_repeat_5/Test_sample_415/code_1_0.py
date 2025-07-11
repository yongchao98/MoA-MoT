import numpy as np
from sklearn.linear_model import LinearRegression

# Step 1: Define the collected data
# Independent Variable 1: Methylation (%) at TSG Promoter
methylation = np.array([10, 15, 20, 25, 30, 35, 40])

# Independent Variable 2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
histone = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent Variable: TSG Expression (RNA-Seq Reads)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Step 2: Prepare the data for the model
# The independent variables (features) are combined into a single matrix X.
X = np.column_stack((methylation, histone))
y = tsg_expression

# Step 3: Create and fit the multiple linear regression model
model = LinearRegression()
model.fit(X, y)

# Step 4: Extract the regression coefficients
beta_0 = model.intercept_      # Intercept
beta_1, beta_2 = model.coef_   # Coefficients for methylation and histone modification

# Step 5: Print the results
print(f"Regression Coefficients:")
print(f"Intercept (β0): {beta_0:.1f}")
print(f"Coefficient for Methylation (β1): {beta_1:.1f}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2:.1f}")

print("\nFinal Regression Equation:")
print(f"TSG Expression = {beta_0:.1f} + ({beta_1:.1f} * Methylation) + ({beta_2:.1f} * Histone H3K9 Trimethylation)")
