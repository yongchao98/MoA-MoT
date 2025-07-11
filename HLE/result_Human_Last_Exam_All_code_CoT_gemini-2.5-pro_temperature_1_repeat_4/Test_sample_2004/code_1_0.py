import numpy as np
from scipy import stats

# The 25 observations of x and y
x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

# The chosen model is of the form y = a*sqrt(x) + b.
# We can find the parameters 'a' and 'b' by performing a linear regression
# on y against the square root of x.
x_transformed = np.sqrt(x_data)

# Perform the linear regression to find the slope (a) and intercept (b)
slope, intercept, _, _, _ = stats.linregress(x_transformed, y_data)

# Format the parameter estimates to 3 significant digits
param_a = f"{slope:.3g}"
param_b_val = abs(intercept)
param_b_formatted = f"{param_b_val:.3g}"

# Determine the sign for the equation string
sign = "+" if intercept >= 0 else "-"

# Print the final model equation, which includes the parameter estimates
print(f"The maximally parsimonious model is:")
print(f"y = {param_a} * sqrt(x) {sign} {param_b_formatted}")