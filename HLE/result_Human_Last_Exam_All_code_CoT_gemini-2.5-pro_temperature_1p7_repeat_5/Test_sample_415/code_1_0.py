import numpy as np
from sklearn.linear_model import LinearRegression

#
# Plan:
# 1. Define the dataset provided by the user using NumPy arrays.
#    - X1: Methylation (%)
#    - X2: Histone H3K9 Trimethylation (ChIP-Seq Reads)
#    - Y: TSG Expression (RNA-Seq Reads)
# 2. Structure the data into a feature matrix (X) for the independent variables
#    and a target vector (y) for the dependent variable.
# 3. Create an instance of the `LinearRegression` model from the scikit-learn library.
# 4. Fit the linear regression model to the data to compute the coefficients.
# 5. Extract the intercept (β0) and the slope coefficients (β1, β2) from the fitted model.
# 6. Print the coefficients and the final regression equation.
#

# Step 1: Define the dataset
methylation = np.array([10, 15, 20, 25, 30, 35, 40])
histone_trimethylation = np.array([300, 275, 250, 225, 200, 175, 150])
tsg_expression = np.array([500, 450, 400, 350, 300, 250, 200])

# Step 2: Structure the data
# The independent variables (features) are combined into a single matrix.
X = np.column_stack((methylation, histone_trimethylation))
# The dependent variable (target).
y = tsg_expression

# Step 3: Create a LinearRegression model instance
model = LinearRegression()

# Step 4: Fit the model to the data
model.fit(X, y)

# Step 5: Extract the coefficients
beta_0 = model.intercept_
beta_1, beta_2 = model.coef_

# Step 6: Print the results
print("The calculated regression coefficients are:")
print(f"Intercept (β0): {beta_0}")
print(f"Coefficient for Methylation (β1): {beta_1}")
print(f"Coefficient for Histone H3K9 Trimethylation (β2): {beta_2}")

print("\nThe final regression equation is:")
# Output each number in the final equation as requested.
print(f"TSG Expression = {beta_0:.4f} + ({beta_1:.4f} * Methylation) + ({beta_2:.4f} * Histone H3K9 Trimethylation)")

# The final answer in the requested format
final_answer = (beta_0, beta_1, beta_2)
print(f"<<<{final_answer}>>>")