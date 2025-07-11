import numpy as np
from scipy import stats

def round_to_sf(num, sf):
    """
    Rounds a number to a specified number of significant figures.
    """
    if num == 0:
        return 0
    # Calculate the position of the most significant digit
    pos = -int(np.floor(np.log10(abs(num))))
    return round(num, pos + sf - 1)

# The 25 observations of x and y
x_data = np.array([
    5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 
    3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45
])
y_data = np.array([
    1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 
    0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 
    2.303, 2.31, 1.194, 1.144, 0.123
])

# Based on analysis, the maximally parsimonious model is y = a * sqrt(x) + b.
# This is equivalent to a linear regression of y vs. sqrt(x).
x_transformed = np.sqrt(x_data)

# Perform linear regression on the transformed data to find parameters 'a' and 'b'.
# slope will be 'a', intercept will be 'b'.
slope, intercept, r_value, p_value, std_err = stats.linregress(x_transformed, y_data)

# Round the parameters to 3 significant digits
a = round_to_sf(slope, 3)
b = round_to_sf(intercept, 3)

# Print the results, including each number in the final equation.
print("The maximally parsimonious model has the form: y = a * sqrt(x) + b")
print(f"The estimated parameter 'a' is: {a}")
print(f"The estimated parameter 'b' is: {b}")

# Determine the sign for the equation string
sign = "+" if b >= 0 else "-"
abs_b = abs(b)

print("\nThe final equation is:")
print(f"y = {a} * sqrt(x) {sign} {abs_b}")

# The final answer content as a single string
final_equation_str = f"y = {a} * sqrt(x) {sign} {abs_b}"