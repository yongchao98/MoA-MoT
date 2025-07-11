import math

# Step 1 & 2: Analyze the properties to determine the homeomorphism class of X.
# The problem states that X is a locally compact metric space and is a
# one-to-one continuous image of the real line R.
# Let f: R -> X be the one-to-one continuous map (a continuous bijection).
#
# A space X is Hausdorff if for any two distinct points, there exist disjoint
# open neighborhoods for them. Metric spaces are always Hausdorff.
#
# So, we have a continuous bijection f between two locally compact Hausdorff spaces
# (R and X). A fundamental theorem in topology states that such a map is a
# homeomorphism.
#
# This implies that X must be topologically equivalent (homeomorphic) to R.

# Step 3: Verify that R satisfies all conditions.
# - R is a metric space (with the standard Euclidean distance).
# - R is locally compact.
# - R is a one-to-one continuous image of itself (via the identity map).
# - We must check the final property: For any two distinct points x, y in R,
#   there exists a closed connected set K such that x is in the interior of K
#   and K is a subset of R \ {y}.
#
# Let x, y be in R, with x != y.
# Let epsilon = abs(x - y) / 2.
# Let K be the closed interval [x - epsilon, x + epsilon].
# K is a closed and connected subset of R.
# The interior of K is the open interval (x - epsilon, x + epsilon), which contains x.
# The distance from any point in K to x is at most epsilon, while the distance from y to x is 2*epsilon.
# Thus, y is not in K. So, K is a subset of R \ {y}.
# All conditions are satisfied by R.

# Step 4: Conclude the number of homeomorphism classes.
# Since any space X satisfying the properties must be homeomorphic to R, and R
# itself is such a space, there is only one possible homeomorphism class.
number_of_classes = 1

# There is no equation in this problem, just a final number to find.
# We print the final answer.
print(number_of_classes)