import numpy as np
from sklearn.linear_model import LinearRegression

# Step 1: Define the dataset
data = np.array([
    [-1, -1, -1, 34.3],
    [ 1, -1, -1, 94.6],
    [-1,  1, -1, 52.5],
    [ 1,  1, -1, 75.1],
    [-1, -1,  1, 28.4],
    [ 1, -1,  1, 48.3],
    [-1,  1,  1,  7.0],
    [ 1,  1,  1, 70.8]
])

# Extract independent (x) and dependent (y) variables
x_vars = data[:, :3]
y_vals = data[:, 3]

# Isolate individual x variables for clarity
x1 = x_vars[:, 0]
x2 = x_vars[:, 1]
x3 = x_vars[:, 2]

# Step 2: Construct the design matrix for the parsimonious model
# Based on preliminary analysis, we select the terms:
# intercept, x1, x3, x2*x3, and x1*x2*x3
X_parsimonious = np.c_[
    np.ones(x_vars.shape[0]),  # Intercept
    x1,
    x3,
    x2 * x3,
    x1 * x2 * x3
]
term_labels = ["", " * x1", " * x3", " * x2*x3", " * x1*x2*x3"]

# Step 3: Fit the linear regression model
model = LinearRegression(fit_intercept=False) # We included the intercept in our matrix
model.fit(X_parsimonious, y_vals)

# Step 4: Get the estimated coefficients and round them
coeffs = [round(c, 1) for c in model.coef_]

# Step 5: Construct and print the final model equation
equation_parts = []
# Handle the intercept first
equation_parts.append(str(coeffs[0]))

# Handle the other terms
for i in range(1, len(coeffs)):
    coeff = coeffs[i]
    label = term_labels[i]
    if coeff != 0:
        sign = "-" if coeff < 0 else "+"
        equation_parts.append(f"{sign} {abs(coeff)}{label}")

final_equation = "y = " + " ".join(equation_parts)
print(final_equation)
print("\nParameter estimates:")
print(f"Intercept: {coeffs[0]}")
print(f"Coefficient for x1: {coeffs[1]}")
print(f"Coefficient for x3: {coeffs[2]}")
print(f"Coefficient for x2*x3: {coeffs[3]}")
print(f"Coefficient for x1*x2*x3: {coeffs[4]}")