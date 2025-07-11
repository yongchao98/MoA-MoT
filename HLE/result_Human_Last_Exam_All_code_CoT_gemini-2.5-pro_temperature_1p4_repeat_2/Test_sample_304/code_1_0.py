# H_max is the maximal entropy we are looking for.
# From the derivation, we have the inequality:
# 2 * H_max = H(y, s1) + H(x, s2)
# Using subadditivity H(A,B) <= H(A)+H(B) and the problem constraints:
# 2 * H_max <= (H(y) + H(s1)) + (H(x) + H(s2))
# We are given H(v) <= 1 for all variables v.
H_y_bound = 1
H_s1_bound = 1
H_x_bound = 1
H_s2_bound = 1

# The upper bound for 2 * H_max is the sum of the individual bounds.
upper_bound_sum = H_y_bound + H_s1_bound + H_x_bound + H_s2_bound

# We print the final equation with the numerical values.
print(f"2 * H_max <= (H(y) + H(s1)) + (H(x) + H(s2)) <= {H_y_bound} + {H_s1_bound} + {H_x_bound} + {H_s2_bound} = {upper_bound_sum}")

# From this, we can find the upper bound for H_max.
H_max_upper_bound = upper_bound_sum / 2
print(f"This implies: H_max <= {int(H_max_upper_bound)}")

# As a valid construction achieving this bound exists, the maximal entropy is this upper bound.
print("The maximal entropy is:")
print(int(H_max_upper_bound))