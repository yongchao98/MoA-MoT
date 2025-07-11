# This script calculates the lower and upper bounds for t based on the problem's requirement.

# Let P be the interval of possible sums, P = [-2, 2t].
# The problem implies that for any x in P, 1/x must also be in P.

# This leads to two main constraints on t.
# First, to avoid division by zero, the sum must never be 0, which means t < 0.
# Second, the interval P must be closed under the reciprocal operation.

# This second constraint gives two inequalities:
# 1. 1/(2*t) >= -2
# 2. -0.5 <= 2*t

# Let's solve these inequalities.

# From inequality 1: 1/(2t) >= -2
# Since t < 0, 2t is negative. Multiplying reverses the inequality.
# 1 <= -4t  =>  t <= -1/4
upper_bound_val = -1.0 / 4.0

# From inequality 2: -0.5 <= 2t
# Dividing by 2: -0.25 <= t
lower_bound_val = -0.5 / 2.0

# We found that t must be both less than or equal to -0.25 and greater than or equal to -0.25.
# This means t can only be -0.25.
t = -0.25

print(f"The first inequality is 1/(2t) >= -2, which simplifies to the final equation: t <= {upper_bound_val}")
print(f"The second inequality is -0.5 <= 2t, which simplifies to the final equation: t >= {lower_bound_val}")
print(f"\nCombining these gives a unique value for t: {t}")
print("The lower and upper bounds are therefore identical.")

# Print the final answer with lower and upper bounds separated by a space.
print(f"\n{lower_bound_val} {upper_bound_val}")