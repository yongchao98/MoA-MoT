# This problem analyzes the configuration space X_4 of 4-segment unit-length
# robot arms that form a closed loop. The space is decomposed into a disjoint
# union of connected manifolds. We are asked to find the dimensions of these manifolds.

# As derived from the analysis, the space X_4 can be stratified by the length L
# of the diagonal vector v1 + v2. This leads to three distinct manifolds:
# 1. The generic case (0 < L < 2), which is a 5-dimensional manifold.
# 2. A degenerate case (L = 0), which is a 4-dimensional manifold (S^2 x S^2).
# 3. Another degenerate case (L = 2), which is a 2-dimensional manifold (S^2).

# The dimensions of these manifolds, ordered from largest to smallest, are y1, y2, y3.
y1 = 5
y2 = 4
y3 = 2

# The final answer is the tuple of these dimensions.
# The code below prints the numbers in the required format.
print(f"{y1},{y2},{y3}")