# The user wants to determine the number of distinct homeomorphism classes
# for a topological space X with specific properties related to the long ray R.

# Step 1: Define the problem mathematically.
# Let R be the long ray, R = [0, omega_1).
# Let X be a compact topological space such that:
# (1) R is a dense subspace of X.
# (2) Every bounded continuous function f: R -> R extends to a unique
#     continuous function F: X -> R.

# Step 2: Analyze the properties.
# Property (1) means that X is a compactification of R.
# Property (2) establishes a special relationship between the function spaces C(X)
# (continuous functions on X) and C_b(R) (bounded continuous functions on R).
#
# Because X is compact, any continuous function F: X -> R is automatically bounded.
# The extension property means that the restriction map, which takes a function F in C(X)
# and restricts its domain to R, is surjective (every f in C_b(R) is an extension of some F).
# The uniqueness of the extension means the restriction map is injective (if two functions
# in C(X) agree on the dense subset R, they must be the same function).
#
# Therefore, property (2) is equivalent to saying that the restriction map
# from C(X) to C_b(R) is an algebra isomorphism: C(X) isomorphic to C_b(R).

# Step 3: Apply the theory of Stone-Čech compactification.
# For any Tychonoff space Y (the long ray R is such a space), its Stone-Čech
# compactification, beta Y, is a compact Hausdorff space that is uniquely defined
# (up to homeomorphism) by the property that every bounded continuous function
# from Y to R extends uniquely to beta Y.
#
# The properties of X given in the problem are precisely the defining
# properties of the Stone-Čech compactification of R.
# Thus, X must be homeomorphic to beta R.

# Step 4: Use the uniqueness property.
# The Stone-Čech compactification is unique up to a homeomorphism that fixes the
# original space. This means that any two spaces, say X_1 and X_2, that both
# satisfy the given conditions must both be homeomorphic to beta R.
# By transitivity of homeomorphism, X_1 and X_2 must be homeomorphic to each other.

# Step 5: Conclude the number of homeomorphism classes.
# Since all spaces satisfying the given conditions belong to the same
# homeomorphism class (the one containing beta R), there is only one such class.
# It is a known theorem in topology that for R = [0, omega_1), its Stone-Čech
# compactification beta R is in fact homeomorphic to the long line [0, omega_1],
# which is the one-point compactification of R.

# Step 6: Final Answer Calculation
# Based on the reasoning, the number of distinct homeomorphism classes is 1.
number_of_classes = 1

# Print the final result.
# The problem asks for the number of distinct homeomorphism classes.
print("The number of distinct homeomorphism classes is:")
print(number_of_classes)
