# The problem asks for the smallest possible cardinality of the set of non-block points
# of an aposyndetic continuum X.
#
# Step 1: We show that for any aposyndetic continuum X, the set of non-block points is X itself.
# This relies on the theorem that an aposyndetic continuum has no cut points, which implies
# that for any p in X, the set X \ {p} is connected. Using the aposyndetic property,
# this can be strengthened to show that X \ {p} is continuum-connected. A point p is
# a non-block point if and only if X \ {p} is continuum-connected.
# Thus, every point in X is a non-block point.
#
# Step 2: The problem is now to find the minimum possible cardinality of an aposyndetic continuum.
#
# Step 3: Consider the degenerate continuum X = {a}, a single-point space.
# - It is a continuum (compact, connected, Hausdorff).
# - It is aposyndetic because the condition is vacuously true (there are no two distinct points).
#
# Step 4: A single-point space is an aposyndetic continuum. Its cardinality is 1.
# Any continuum with more than one point must be infinite. Therefore, the minimum cardinality is 1.
# The set of non-block points for X = {a} is {a}, which has cardinality 1.

smallest_cardinality = 1
print(smallest_cardinality)