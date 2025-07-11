import numpy as np

# The limit of M(theta) as theta -> 0 is given by arccos(2*sqrt(6)/5).
# We calculate this value.

num_2 = 2
num_6 = 6
num_5 = 5

# Calculate the value of the cosine of the angle.
cos_M = num_2 * np.sqrt(num_6) / num_5

# Calculate the angle M in radians.
M = np.arccos(cos_M)

print("The limit of M(theta) as theta approaches 0 is given by the expression arccos(2 * sqrt(6) / 5).")
print(f"The number in the expression are: {num_2}, {num_6}, and {num_5}.")
print(f"The numerical value of this limit is: {M}")

# We can also express the answer as arcsin(1/5)
M_arcsin = np.arcsin(1/5)
# print(f"This is equivalent to arcsin(1/5), which is: {M_arcsin}")
