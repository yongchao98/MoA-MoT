import numpy as np
from scipy import stats

# The 25 observations of x and y
x_obs = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
y_obs = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

# The chosen model is of the form y = a*sqrt(x) + b.
# We can find the parameters 'a' and 'b' using a linear regression on y vs. sqrt(x).
x_transformed = np.sqrt(x_obs)

# Perform the linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x_transformed, y_obs)

# The slope of this regression is our 'a' and the intercept is our 'b'.
a = slope
b = intercept

# Format the parameters to 3 significant digits.
# The formatg function handles scientific notation and rounding for us.
a_str = "{:.3g}".format(a)
b_str = "{:.3g}".format(b)

# Print the final model equation.
# The sign of b is handled by the f-string formatting.
print(f"The optimal model is y = {a_str} * sqrt(x) {b:+.3g}")

# As per the instruction "Remember in the final code you still need to output each number in the final equation!"
# The line above prints the full equation. Let's break down the numbers themselves.
# print(f"Parameter 'a': {a_str}")
# print(f"Parameter 'b': {b_str}")