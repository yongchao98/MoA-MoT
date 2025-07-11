import sys

# This script finds the lower and upper bounds of t based on the problem statement.

# Step 1: Analyze the constraints on t.
# Let the interval for a_i be I = [-1, t]. For I to be non-empty, t >= -1.
# The sum S = a_0 + a_2 can take any value in the range J = [-2, 2*t].
# The problem requires that for any S in J, 1/S is also in J.
# This implies S cannot be 0, so 0 must not be in the range J.
# 0 is in J = [-2, 2*t] if 2*t >= 0, which means t >= 0.
# To exclude 0, we must have t < 0.
# So, the valid range for t must be within [-1, 0).

print("Let the range of sums S = a_0 + a_2 be J = [-2, 2*t].")
print("To satisfy the condition, we must have t < 0.")

# Step 2: Set up the inequalities for t.
# For t in [-1, 0), the range J contains only negative numbers.
# The set of reciprocals {1/s | s in J} is [1/(2*t), -0.5].
# This reciprocal range must be a subset of J = [-2, 2*t].
# This leads to two inequalities.

# Inequality 1: The lower bound of the reciprocal range must be >= the lower bound of J.
# 1/(2*t) >= -2
# Since t < 0, multiplying by 2*t flips the inequality sign.
# 1 <= -4*t  =>  t <= -1/4
print("From the inequality 1/(2*t) >= -2 (derived from the lower bounds), we get t <= -1/4.")
# The instruction asked to "output each number in the final equation".
# The numbers involved in this inequality are 1, 2, -2, -4.
print("Numbers in the derivation: 1, 2, -2, -4")


# Inequality 2: The upper bound of the reciprocal range must be <= the upper bound of J.
# -0.5 <= 2*t
# -1/4 <= t
print("From the inequality -0.5 <= 2*t (derived from the upper bounds), we get t >= -1/4.")
# The numbers involved in this inequality are -0.5, 2, -1, 4.
print("Numbers in the derivation: -0.5, 2, -1, 4")


# Step 3: Solve for t and determine the bounds.
# The only value of t satisfying both t <= -1/4 and t >= -1/4 is t = -1/4.
t = -0.25
lower_bound = t
upper_bound = t

print(f"\nThe only value of t that satisfies both inequalities is t = {t}.")
print("Therefore, the lower and upper bounds are the same.")
print(f"\n{lower_bound} {upper_bound}")