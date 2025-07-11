import numpy as np
from sklearn.linear_model import LinearRegression

# Step 1: I am setting up the data from the table into NumPy arrays.
# Independent variables (predictors)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
h3k9_trimethylation = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (response)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Step 2: I will combine the two predictor variables into a single feature matrix.
# scikit-learn expects the features in a (n_samples, n_features) format.
X = np.column_stack((methylation, h3k9_trimethylation))
y = tsg_expression

# Step 3: I will create an instance of the Linear Regression model and fit it to the data.
model = LinearRegression()
model.fit(X, y)

# Step 4: I will extract the calculated coefficients from the trained model.
# The intercept (β0)
beta_0 = model.intercept_
# The slope coefficients (β1 for methylation, β2 for H3K9 trimethylation)
beta_1, beta_2 = model.coef_

# The calculated coefficients may have minor floating-point inaccuracies, so I will round them for a clean presentation.
b0_final = round(beta_0)
b1_final = round(beta_1)
b2_final = round(beta_2)

# Step 5: I will print the results in a clear and understandable format.
print("Multiple Linear Regression Results:")
print(f"Intercept (β0): {b0_final}")
print(f"Coefficient for Methylation (β1): {b1_final}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {b2_final}")

print("\nThe final regression equation is:")
print(f"TSG Expression = {b0_final} + ({b1_final} * Methylation) + ({b2_final} * Histone H3K9 Trimethylation)")
