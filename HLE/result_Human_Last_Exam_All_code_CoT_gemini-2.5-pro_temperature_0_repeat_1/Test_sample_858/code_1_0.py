# Step 1: Relate non-block points to a simpler concept.
# A point p in a continuum X is a non-block point if X \ {p} contains a dense continuum-connected subset.
# A point p in a continuum X is a non-cut-point if X \ {p} is connected.
#
# For an aposyndetic continuum, we can prove that the set of non-block points is exactly the set of non-cut-points.
# - If p is a non-cut-point of an aposyndetic continuum, X \ {p} is connected. It can be shown that this implies X \ {p} is also continuum-connected. A continuum-connected space is its own dense continuum-connected subset. Thus, p is a non-block point.
# - Conversely, if p is a cut point, X \ {p} is disconnected. Any continuum-connected subset of X \ {p} must lie entirely within one of the connected components of X \ {p}. Such a subset cannot be dense in the entire space X \ {p}. Therefore, a cut point is always a block point.
#
# So, the problem is equivalent to finding the smallest possible cardinality of the set of non-cut-points of an aposyndetic continuum.

# Step 2: Analyze cases for the continuum X to find the minimum number of non-cut-points.
# We consider two cases: degenerate (a single point) and non-degenerate (more than one point) continua.

# Step 3: Case 1: The degenerate continuum.
# Let X be a continuum consisting of a single point, X = {p}.
# - Is X a continuum? Yes, it is compact, connected, and Hausdorff.
# - Is X aposyndetic? The condition for aposyndesy is "for every two distinct points x, y...". Since there are no two distinct points in X, the condition is vacuously true. So, X is aposyndetic.
# - How many non-cut-points does X have? A point q is a non-cut-point if X \ {q} is connected. For q=p, X \ {p} is the empty set. The empty set is connected.
# - Therefore, p is a non-cut-point. The set of non-cut-points is {p}.
# - The cardinality of the set of non-cut-points (and thus non-block points) is 1.

# Step 4: Case 2: The non-degenerate continuum.
# Let X be a continuum with more than one point.
# - A classical theorem by R. L. Moore states that any non-degenerate continuum has at least two non-cut-points.
# - Therefore, any non-degenerate aposyndetic continuum must have at least 2 non-block points.
# - An example that achieves this minimum is the closed interval X = [0, 1]. It is an aposyndetic continuum. Its cut points are all points in (0, 1). Its non-cut-points are the endpoints {0, 1}.
# - The cardinality of the set of non-block points for X = [0, 1] is 2.

# Step 5: Conclusion.
# Comparing the two cases, the minimum cardinality for the set of non-block points is 1 (from the degenerate case). The problem statement does not exclude the degenerate continuum.
# The smallest possible cardinality is 1.

# The final equation is simply the result of our logical deduction.
final_answer = 1
print(f"The smallest possible cardinality of the set of non-block points is: {final_answer}")