import numpy as np

# Based on the derivation, the only value for t that satisfies the condition is t = -1/4.
# Let S be the interval of possible sums [a_0 + a_2], which is [-2, 2t].
# The condition implies that for any s in S, 1/s must also be in S.
# First, to avoid division by zero, 0 cannot be in S, which means 2t < 0, so t < 0.
# Let I = [-2, 2t]. The image of I under f(x)=1/x is [1/(2t), -1/2].
# We require [1/(2t), -1/2] to be a subset of [-2, 2t].
# This gives two inequalities:
# 1) -2 <= 1/(2t)  => -4t >= 1 (since t<0) => t <= -1/4
# 2) -1/2 <= 2t => t >= -1/4
# Combining t <= -1/4 and t >= -1/4, we get t = -1/4.

lower_bound = -0.25
upper_bound = -0.25

print(f"{lower_bound} {upper_bound}")