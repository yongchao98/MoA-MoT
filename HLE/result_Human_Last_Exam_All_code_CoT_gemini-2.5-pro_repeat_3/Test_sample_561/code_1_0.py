import math

# N is the number of self-similar copies at each iteration.
# A one-octave keyboard has 5 black keys.
N = 5

# s is the scaling factor.
# The original keyboard (3x1) must be scaled to fit onto a black key (3/14 x 9/14).
# The scaling factor s must satisfy:
# 3*s <= 3/14  => s <= 1/14
# 1*s <= 9/14  => s <= 9/14
# The more restrictive (smaller) scaling factor is s = 1/14.
s = 1/14

# The Minkowski-Bouligand dimension D is log(N) / log(1/s)
dimension = math.log(N) / math.log(1/s)

# Print the final equation and the result
print(f"The Minkowski-Bouligand dimension D is calculated as log(N) / log(1/s).")
print(f"D = log({N}) / log(1/({s:.4f}))")
print(f"D = log({N}) / log({1/s})")
print(f"D = {dimension}")