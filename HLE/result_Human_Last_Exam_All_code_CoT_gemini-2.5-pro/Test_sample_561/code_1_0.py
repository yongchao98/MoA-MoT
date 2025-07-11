import math

# N is the number of self-similar copies. In one octave, there are 5 black keys.
N = 5

# r is the scaling factor.
# The original keyboard is 3 units wide. It has 7 white keys.
# Width of a white key = 3 / 7
# Width of a black key = (3 / 7) / 2 = 3 / 14
# The new keyboard has width 3*r, which must be <= 3/14. So r <= 1/14.
#
# The original keyboard is 1 unit high.
# Height of a black key = 9 / 14
# The new keyboard has height 1*r, which must be <= 9/14. So r <= 9/14.
#
# To fit, we must take the more restrictive (smaller) scaling factor.
# r = 1/14
r = 1/14

# The Minkowski–Bouligand dimension D is log(N) / log(1/r)
# D = log(5) / log(1/(1/14)) = log(5) / log(14)
s = 1/r # s is the magnification factor
dimension = math.log(N) / math.log(s)

# Print the final equation with the numbers used
print(f"The dimension D is calculated as log(N) / log(1/r)")
print(f"D = log({N}) / log(1/({r:.4f}))")
print(f"D = log({N}) / log({s})")
print(f"D = {dimension}")
