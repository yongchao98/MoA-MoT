# The problem is to determine the number of distinct homeomorphism classes for a space X
# with specific properties.
#
# Let X be a compact, connected metric space.
# Let C_n(X) be the configuration space of n distinct points from X.
# We are given that for some n >= 2, C_n(X) is disconnected.
#
# A key theorem in topology states that for a non-degenerate compact connected metric space X,
# the configuration space C_n(X) is disconnected for some n >= 2 if and only if X is
# homeomorphic to the closed interval [0, 1].
#
# This means that any space X satisfying the given conditions must be an "arc".
# All arcs are homeomorphic to each other. For example, the interval [0, 1] is
# homeomorphic to [a, b] for any a < b, and also to any simple, non-self-intersecting
# curve with two endpoints.
#
# Since all such spaces X belong to the same homeomorphism class (the class of the
# interval [0, 1]), there is only one such class.

# The number of distinct homeomorphism classes.
number_of_homeomorphism_classes = 1

# Print the final result.
print(number_of_homeomorphism_classes)