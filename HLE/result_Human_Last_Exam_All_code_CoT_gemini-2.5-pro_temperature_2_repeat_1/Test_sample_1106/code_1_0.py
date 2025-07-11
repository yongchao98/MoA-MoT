# Number of points
n = 15
# Dimension of the ambient Euclidean space
d = 8

# The problem is to find the minimized maximum number of points in any closed
# hemisphere for n points on a (d-1)-sphere.
# For the special case where n = 2d - 1 and d is a power of two,
# this value is known to be exactly d.

# The given values n=15 and d=8 satisfy these conditions.
# The result k is therefore d, which can be derived from n as follows:
k = (n + 1) / 2

# We output the final equation, showing that the result is d and
# how it can be calculated from n.
# The equation format is "d = (n + 1) / 2 = k"
print(f"{d} = ({n} + 1) / 2 = {int(k)}")