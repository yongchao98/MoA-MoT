import numpy as np
from scipy import stats

def round_to_sf(num, sf):
    """
    Rounds a number to a specified number of significant figures.
    """
    if num == 0:
        return 0
    # Calculate the position of the first significant digit
    pos = -int(np.floor(np.log10(abs(num))))
    # Round the number to the correct number of decimal places
    return round(num, pos + sf - 1)

# The 25 observations of x and y
x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

# Based on analysis, the chosen model is y = a * sqrt(x) + b.
# We linearize this by letting z = sqrt(x), so the model becomes y = a*z + b.
z_data = np.sqrt(x_data)

# Perform linear regression to find the best-fit parameters a (slope) and b (intercept).
lin_reg_result = stats.linregress(z_data, y_data)

# Extract the slope (a) and intercept (b)
a = lin_reg_result.slope
b = lin_reg_result.intercept

# Round the parameters to 3 significant digits
a_rounded = round_to_sf(a, 3)
b_rounded = round_to_sf(b, 3)

# Prepare the final equation string for printing
if b_rounded >= 0:
    sign = "+"
    b_val = b_rounded
else:
    sign = "-"
    b_val = -b_rounded

# Print the final equation with the estimated parameters
print("The maximally parsimonious model for predicting y from x is:")
print(f"y = {a_rounded} * sqrt(x) {sign} {b_val}")