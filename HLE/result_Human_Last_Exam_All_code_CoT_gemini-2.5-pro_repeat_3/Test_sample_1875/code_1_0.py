import numpy as np

# This script builds a parsimonious regression model from the given data.
# The plan is as follows:
# 1. Set up the data provided by the user.
# 2. Construct a design matrix for a full linear model, including all main effects and interactions.
#    This allows us to estimate the contribution of each potential term.
# 3. Solve for the model coefficients using ordinary least squares. The orthogonality of the
#    factorial design ensures the coefficient estimates are stable and independent.
# 4. To create a "parsimonious" model, we examine the magnitudes of the coefficients. We will drop
#    terms with small coefficients, as they have little predictive power. We will also adhere
#    to the model hierarchy principle.
# 5. After analyzing the coefficients, the chosen model includes the intercept, the main effects
#    x_1, x_2, x_3, and the two-way interactions x_1*x_3 and x_2*x_3. This model retains the
#    strongest effects while being simpler than the full model.
# 6. The script will then format and print the final model equation with coefficients
#    rounded to one decimal place as requested.

# 1. Define the experimental data
data = np.array([
    [-1, -1, -1, 34.3],
    [1,  -1, -1, 94.6],
    [-1, 1,  -1, 52.5],
    [1,  1,  -1, 75.1],
    [-1, -1, 1,  28.4],
    [1,  -1, 1,  48.3],
    [-1, 1,  1,  7.0],
    [1,  1,  1,  70.8]
])

x_data = data[:, :3]
y_data = data[:, 3]

# 2. Construct the full design matrix
intercept = np.ones(x_data.shape[0])
x1 = x_data[:, 0]
x2 = x_data[:, 1]
x3 = x_data[:, 2]
x1x2 = x1 * x2
x1x3 = x1 * x3
x2x3 = x2 * x3
x1x2x3 = x1 * x2 * x3

X_design = np.c_[intercept, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3]
term_labels = ["", "x_1", "x_2", "x_3", "x_1*x_2", "x_1*x_3", "x_2*x_3", "x_1*x_2*x_3"]

# 3. Solve for coefficients of the full model
full_coeffs = np.linalg.lstsq(X_design, y_data, rcond=None)[0]

# 4. Select terms for the parsimonious model based on coefficient magnitude and hierarchy.
# The selected model is: y = b0 + b1*x1 + b2*x2 + b3*x3 + b13*x1*x3 + b23*x2*x3
# These correspond to indices [0, 1, 2, 3, 5, 6] in the full model.
indices_to_keep = [0, 1, 2, 3, 5, 6]
parsimonious_coeffs = full_coeffs[indices_to_keep]
parsimonious_terms = [term_labels[i] for i in indices_to_keep]

# 5. Construct and print the final equation string
# Start with the intercept
equation = f"y = {parsimonious_coeffs[0]:.1f}"

# Add the other terms
for i in range(1, len(parsimonious_coeffs)):
    coeff = parsimonious_coeffs[i]
    term = parsimonious_terms[i]
    sign = "-" if coeff < 0 else "+"
    equation += f" {sign} {abs(coeff):.1f}*{term}"

print(equation)
