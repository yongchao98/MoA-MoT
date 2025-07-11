import numpy as np

def round_to_sf(num, sf):
    """
    Rounds a number to a specified number of significant figures.
    """
    if num == 0:
        return 0
    # Calculate the power of 10 to round to
    round_pos = sf - int(np.floor(np.log10(abs(num)))) - 1
    return round(num, round_pos)

# 1. Define the dataset
x_obs = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
y_obs = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

# 2. Choose the best model y = c + a*sqrt(x) and perform a linear fit
# We transform x to sqrt(x) and then fit a line between y and sqrt(x).
sqrt_x = np.sqrt(x_obs)

# np.polyfit with degree 1 performs a linear regression.
# It returns coefficients [slope, intercept] which correspond to [a, c] in our model.
coeffs = np.polyfit(sqrt_x, y_obs, 1)
a = coeffs[0]
c = coeffs[1]

# 3. Round the parameters to 3 significant figures
a_rounded = round_to_sf(a, 3)
c_rounded = round_to_sf(c, 3)

# 4. Print the final equation with the calculated parameters
# The numbers in the final equation are a_rounded and c_rounded.
# We format the string to show the negative sign correctly.
print("The determined maximally parsimonious model is of the form: y = a * sqrt(x) + c")
print("\nThe estimated parameters are:")
print(f"a = {a_rounded}")
print(f"c = {c_rounded}")
print("\nThe final equation is:")
print(f"y = {a_rounded} * sqrt(x) - {abs(c_rounded)}")