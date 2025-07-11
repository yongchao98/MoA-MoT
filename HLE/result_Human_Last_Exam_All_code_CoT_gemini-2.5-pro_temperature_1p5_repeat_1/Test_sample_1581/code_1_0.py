# The problem asks for the number of distinct homeomorphism classes of spaces with a certain property.
# Let X be a compact connected metric space.
# The property is that for some integer n >= 2, the configuration space C_n(X) of n distinct points in X is disconnected.

# A classical theorem in continuum theory provides a full characterization for such spaces.
# Theorem: A compact connected metric space X is homeomorphic to the closed interval [0, 1]
# if and only if C_n(X) is disconnected for some n >= 2.

# This means that any space X that satisfies the conditions given in the problem
# must be homeomorphic to the interval [0, 1].

# The question is to find the number of distinct homeomorphism classes.
# All spaces that are homeomorphic to each other belong to the same homeomorphism class.
# Since any such space X must be homeomorphic to [0, 1], all such spaces fall into
# a single homeomorphism class (the one containing [0, 1]).

# Therefore, there is only one such homeomorphism class.

# Let N be the number of distinct homeomorphism classes.
N = 1

# The problem asks to output the final equation.
print(f"Number of distinct homeomorphism classes = {N}")