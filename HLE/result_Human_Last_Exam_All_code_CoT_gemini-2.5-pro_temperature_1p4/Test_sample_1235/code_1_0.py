import numpy as np

# Define the constants from the equation for the generating amplitude c_1
# The equation is: 2*pi - (c_1**2 / 2) * (e**(4*pi) - 1) = 0
num_2 = 2
num_4 = 4
num_1 = 1

# Numerator of the expression for c_1^2 is 4*pi
numerator = num_4 * np.pi
# Denominator of the expression for c_1^2 is e^(4*pi) - 1
denominator = np.exp(num_4 * np.pi) - num_1

# Calculate c_1^2
c1_squared = numerator / denominator

# c_1 must be positive
c1 = np.sqrt(c1_squared)

print("The equation for the generating amplitude c_1 is:")
print(f"{num_2}*pi - (c_1^2 / {num_2}) * (e^({num_4}*pi) - {num_1}) = 0")
print("\nRearranging for c_1^2 gives:")
print(f"c_1^2 = ({num_4}*pi) / (e^({num_4}*pi) - {num_1})")
print(f"c_1^2 = {numerator} / ({np.exp(num_4 * np.pi)} - {num_1})")
print(f"c_1^2 = {c1_squared}")
print("\nThe value of the first positive root c_1 is:")
print(f"c_1 = {c1}")