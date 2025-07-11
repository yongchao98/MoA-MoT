import numpy as np
from sklearn.linear_model import LinearRegression

# 1. Load the Data
# Independent variables (features)
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_trimethylation = np.array([300, 275, 250, 225, 200, 175, 150])

# Dependent variable (target)
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# 2. Prepare the Data for Regression
# Combine independent variables into a single feature matrix X
# The shape needs to be (n_samples, n_features)
X = np.vstack((methylation, histone_trimethylation)).T
y = tsg_expression

# 3. Perform Multiple Linear Regression
# Create and fit the model
model = LinearRegression()
model.fit(X, y)

# 4. Extract Coefficients
# β0 is the intercept
beta_0 = model.intercept_
# β1 and β2 are the slope coefficients
beta_1, beta_2 = model.coef_

# 5. Display the Results
print("Multiple Linear Regression Coefficients:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")
print("\nFinal Regression Equation:")
# The formatting ensures we show the signs correctly
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * H3K9_Trimethylation)")

# The data in this problem exhibits perfect multicollinearity, meaning the independent
# variables are perfectly correlated. This can lead to unstable coefficient estimates.
# The solution provided is one of many possible mathematical solutions.
# For example, another valid equation is: TSG Expression = 600.0 - 10.0 * Methylation + 0.0 * H3K9_Trimethylation
# The solution from scikit-learn is chosen based on its internal algorithm.

final_answer = f"{beta_0},{beta_1},{beta_2}"
print(f"\n<<<__{final_answer}__>>>")