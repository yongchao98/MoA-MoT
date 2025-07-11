# This problem is a mathematical one concerning the topology of a configuration space.
# The solution involves analyzing the structure of the space X_4 of closed
# 4-segment unit-length arms.
# The space can be stratified based on the length 'r' of the diagonal connecting
# the start of the first segment to the end of the second. This length r can
# range from 0 to 2.
# This stratification leads to a decomposition of X_4 into three disjoint
# connected manifolds.

# 1. The generic case: 0 < r < 2. The dimension of this manifold is 5.
y1 = 5

# 2. The degenerate case: r = 0. This corresponds to configurations of the form
# (v1, -v1, v3, -v3), which is the space S^2 x S^2. The dimension is 4.
y2 = 4

# 3. The degenerate case: r = 2. This corresponds to configurations of the form
# (v1, v1, -v1, -v1), which is the space S^2. The dimension is 2.
y3 = 2

# The problem asks for the tuple of dimensions, ordered from largest to smallest.
# The final answer is the sequence of these dimensions.
print(f"{y1},{y2},{y3}")