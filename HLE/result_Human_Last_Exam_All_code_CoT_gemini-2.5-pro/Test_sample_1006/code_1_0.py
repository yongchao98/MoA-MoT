# The problem asks for the number of distinct homeomorphism classes for a compact topological space X
# with specific properties related to the long ray R.

# The properties are:
# 1. X contains a dense copy of the long ray R. This means X is a compactification of R.
# 2. Every bounded continuous function f: R -> R extends to a unique continuous function on X.

# These two properties are the defining characteristics of the Stone-Čech compactification of R, denoted as βR.
# A key theorem in topology states that the Stone-Čech compactification of a space is unique up to homeomorphism.

# This means that any space X satisfying the given conditions must be homeomorphic to βR.
# Since all such spaces are homeomorphic to each other, they all belong to the same homeomorphism class.

# Therefore, there is only one such distinct homeomorphism class.

number_of_classes = 1

# The problem asks to output the number in the final equation.
# The final equation is simply: Number of classes = 1.
print("The number of distinct homeomorphism classes is:")
print(number_of_classes)