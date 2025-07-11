import math

# Based on the derived inequalities:
# n_r < 3
# floor(n_y / 2) <= n_r
# floor(n_g / 2) <= n_y

# To maximize the total n = n_r + n_g + n_y, we select the largest possible values
# for each variable that satisfy the constraints, starting with n_r.
# The largest integer value for n_r less than 3 is 2.
n_r = 2

# Now, we find the maximum n_y satisfying floor(n_y / 2) <= 2.
# This inequality is equivalent to n_y / 2 < 3, which means n_y < 6.
# The largest integer n_y is 5.
n_y = 5

# Finally, we find the maximum n_g satisfying floor(n_g / 2) <= 5.
# This inequality is equivalent to n_g / 2 < 6, which means n_g < 12.
# The largest integer n_g is 11.
n_g = 11

# The maximum counts for the three colors are a permutation of (2, 5, 11).
# We sort them for a canonical representation.
counts = sorted((n_r, n_y, n_g))
c1, c2, c3 = counts[0], counts[1], counts[2]

# The maximum value of n is the sum of these counts.
max_n = c1 + c2 + c3

print(f"The maximum number of points n is found by summing the largest possible counts for each color, which are {c1}, {c2}, and {c3}.")
print("The final equation for the maximum value is:")
print(f"{c1} + {c2} + {c3} = {max_n}")
