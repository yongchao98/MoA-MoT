import math

# This script calculates and prints the range for the median m.
# The logic is based on the geometric properties derived in the plan.

# Lower bound for m
lower_bound = 13

# The upper bound calculation
# m < sqrt(144 + (1440/119)^2)
# m < sqrt(12^2 * (1 + (120/119)^2))
# m < 12 * sqrt((119^2 + 120^2)/119^2)
# m < (12/119) * sqrt(14161 + 14400)
# m < (12/119) * sqrt(28561)
# sqrt(28561) is 169
# m < (12 * 169) / 119
upper_bound_num = 12 * 169
upper_bound_den = 119

# Print the final inequality with each number.
print("The range of values for m for which angle A is acute is:")
print(f"{lower_bound} < m < {upper_bound_num}/{upper_bound_den}")